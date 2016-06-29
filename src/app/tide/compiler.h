// Benjamin Diament

// We want to be able to compute a dot product between a theoretical peak set
// and an observed peak set. Under the assumption (borne out by use cases seen
// in many test sets) that the same peptide may be a candidate for many
// different observed spectra, we generate machine code, on the fly, for taking
// the dot product between a fixed theoretical peak set and a variable observed
// peak set. 
//
// That is, the generated code is specific to the theoretical peak set, but not
// to the observed peak set. A program is generated for a given peptide's
// theoretical peak set and remains fixed. This program may be run on many
// different inputs representing many different observed spectra.  The input to
// the generated program consists of the cached vectors for an observed peak set
// (as described in spectrum_preprocess.h).
//
// The generated program has a very simple form: it simply adds together certain
// specific cache entries specified by an index, and subtracts others. The
// correct cache index is given by a specific theoretical peak. The cache
// contents will have been determined during spectrum preprocessing (see
// spectrum_preprocess.h).
//
// The programs for all peptides in the ActivePeptideQueue appear
// consecutively in memory, (but see caveat about FifoAllocator below) and are
// all run, one right after the other, on the same input.  The input is a
// buffer pointing to the cache for the observed spectrum. Each program writes
// an entry to an output buffer with a pair (score, counter).  A counter is
// initialized with the number of candidate peptides to score.
//
// The generated code is roughly equivalent to the following C code (the array
// indices are for illustration, the actual values are determined by the code
// generator for the particluar candidate peptide's theoretical peak set):
//
//    *outputbuf++ = cache[72] + cache[195] - cache[546] + ...
//    *outputbuf++ = counter;
//    if (--counter == 0) return;
//    *outputbuf++ = cache[90] + cache[108] + cache[125] + ...
//    *outputbuf++ = counter;
//    if (--counter == 0) return;
//    *outputbuf++ = cache[32]] - cache[42] - cache[321] + ...
//    *outputbuf++ = counter;
//    if (--counter == 0) return;
//    ...
//
// As new peptides enter the ActivePeptideQueue the above program is extended.
// As "older" (lighter) peptides are dequeued and discarded corresponding
// early lines of the generated program are discarded. Each Peptide maintains
// a pointer to begining of its corresponding line of generated code. A
// FifoAllocator manages memory for the generated instructions. When the end
// of a FifoAllocator page is reached, a jump instruction is inserted to
// continue execution at the next FifoAllocator page.
//
// The compiler here generates x86 code, and assumes cache entries are 32-bit
// integers. Each program assumes the following upon entry:
//
//    Register ECX contains the counter i.e. the number of peptides yet to be
//    scored.  When ECX = 0 control returns to the caller.
//
//    Register EDI (RDI on 64-bit machines) points to a buffer where the
//    output of the current program will go (outputbuf in the illustrative C
//    code above).
//
//    Register EDX points to the cached vectors for an observed peak set
//    (cache in the illustrative C code above). 
//
// For illustration the first three lines of the above C code would be compiled as follows:
//
//    mov (%edx+288), %eax // EAX = cache[72] (since 288 == 72 * 4)
//    add (%edx+780), %eax // EAX += cache[195] (since 780 == 195 * 4)
//    sub (%edx+2184), %eax // EAX -= cache[546] (since 2184 == 546 * 4)
//    ...
//    stosl // equivalent to mov %eax, (%edi); inc %edi; (*outputbuf++=%eax)
//    mov %ecx, %eax
//    stosl // equivalent to mov %eax, (%edi); inc %edi; (*outputbuf++=%ecx)
//    loop +1 // equivalent to dec %ecx; if (ecx != 0) skip one instruction
//    ret
//    ... (next program here)

#ifndef COMPILER_H
#define COMPILER_H

#include <stdint.h>

class TheoreticalPeakCompiler {
 public:
  explicit TheoreticalPeakCompiler(FifoAllocator* fifo_alloc) 
    : fifo_alloc_(fifo_alloc), last_alloc_end_(NULL) {
      // fifo_alloc_ will make room for generated programs.
  }

  void* Init(int pos_size, int neg_size) {
    // Init() gets called once per candidate peptide.
    // pos_size is the number of cache entries to be added together, neg_size
    // is the number to be subtracted.
    // add or sub instructions take six bytes. The coda (containing the
    // storage of results etc.) takes seven bytes.
    int total_size = 6*(pos_size + neg_size) + 7;
    pos_ = (unsigned char*) fifo_alloc_->New(total_size); 
    // last_alloc_end points just beyond the last allocated program. By the
    // end of this block, we will ensure that there is enough room for a five
    // byte jump instruction in case current program didn't end up being
    // allocated right after the old.
    if (pos_ == last_alloc_end_) { // Are programs consecutive?
      last_alloc_end_ = pos_ + total_size;
      pos_ -= jmp_size; // Last prog. didn't need room for jump instruction.
    } else { // Progs. not consecutive. Old prog. needs to jump to new prog.
       // Above allocation wasn't enough. Need room for jump.
      fifo_alloc_->Unalloc(total_size);
      total_size += jmp_size;
      pos_ = (unsigned char*) fifo_alloc_->New(total_size);
      if (last_alloc_end_ != NULL)
        AddJump(last_alloc_end_ - jmp_size, pos_);
      last_alloc_end_ = pos_ + total_size;
    }
    first_ = true; // first theoretical peak gets handled a bit differently. 
    return pos_; 
  }

  // The following three functions assume peaks unsorted, and check that each peak
  // doesn't point past end of cache.
    
  void AddPositive(const TheoreticalPeakArr& peaks) {
    // Write an add instruction for each entry in peaks.
    int end = MaxBin::Global().CacheBinEnd() * NUM_PEAK_TYPES;
    for (int i = 0; i < peaks.size(); ++i)
      if (peaks[i].Code() < end)
        AddPositive(peaks[i].Code());
  }

  void AddPositive(const google::protobuf::RepeatedField<int>& peaks) {
    // Write an add instruction for each entry in peaks.
    int end = MaxBin::Global().CacheBinEnd() * NUM_PEAK_TYPES;
    int total = 0;
    google::protobuf::RepeatedField<int>::const_iterator i = peaks.begin();
    for (; i != peaks.end(); ++i) {
      if ((total += *i) >= end)
        break;
      AddPositive(total);
    }
  }

  void AddNegative(const google::protobuf::RepeatedField<int>& peaks) {
    // Write a sub instruction for each entry in peaks.
    int end = MaxBin::Global().CacheBinEnd() * NUM_PEAK_TYPES;
    int total = 0;
    google::protobuf::RepeatedField<int>::const_iterator i = peaks.begin();
    for (; i != peaks.end(); ++i) {
      if ((total += *i) >= end)
        break;
      AddNegative(total);
    }
  }

  void Done() {
    // Write the coda instructions which will store results and update
    // counter. See comments above.
    // Poke machine code into the next 7 bytes at pos_.
    // (int*) pos_ must be a four byte pointer!
    *((int*) pos_) = 0xc889ab; // stosl; mov %ecx, %eax
    *((int*) (pos_+3)) = 0xc301e2ab; // stosl; loop +1; ret
    assert(pos_ + 7 + jmp_size <= last_alloc_end_); // There should still be room
                                                    // for a jump if needed
    last_alloc_end_ = pos_ + 7 + jmp_size;
    fifo_alloc_->Unalloc(last_alloc_end_);
  }

 private:
  // Various x86 instructions we need.
  static const uint16_t add_to_eax_at_edx_plus = 33283;
  static const uint16_t mov_to_eax_at_edx_plus = 33419;
  static const uint16_t sub_from_eax_at_edx_plus = 33323;
  static const unsigned char ret = 195;
  static const unsigned char jmp_relative = 233;

  static const int jmp_size = 5;

  void AddPositive(int peak) {
    if (first_) { // First theoretical peak uses a 'mov' rather than an 'add'
      *((uint16_t*) pos_) = mov_to_eax_at_edx_plus;
    } else {
      *((uint16_t*) pos_) = add_to_eax_at_edx_plus;
    }
    first_ = false;
    pos_ += 2;
    *((int*) pos_) = peak << 2; // Store 4 * the peak position.
    pos_ += 4;
  }

  void AddNegative(int peak) {
    *((uint16_t*) pos_) = sub_from_eax_at_edx_plus;
    pos_ += 2;
    *((int*) pos_) = peak << 2; // Store 4 * the peak position.
    pos_ += 4;
  }

  static void AddJump(unsigned char* pos, unsigned char* whereto) {
    int diff = whereto - (pos + jmp_size);
    *pos++ = jmp_relative;
    *((int*) pos) = diff;
  }

  FifoAllocator* fifo_alloc_;
  unsigned char* last_alloc_end_;
  unsigned char* pos_; // "cursor position" as we write out instructions.
  bool first_;
};

#endif // COMPILER_H
