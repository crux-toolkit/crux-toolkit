// Benjamin Diament
//
// Memory allocator for FIFO usage pattern: client may allocate arbitrarily 
// sized blocks sequentially. A pointer p is returned to each new block.
// Subsequently a call to Release(p) frees all allocations prior to (not including)
// allocation of p. It's not possible to free blocks in any order other than that
// in which they were allocated. This prevents fragmentation but is only appropriate
// for FIFO usage patterns, e.g. data assocatied with a queue.
// 
// A page size, S, is supplied to the FifoAllocator constructor.
// At most 2 * S extra memory will be allocated.
//
// Not thread safe! (TODO 254)
//
// Unalloc() allows you to deallocate the most recently allocated pointer.
//
// Example usage:
// FifoAllocator alloc(1 << 20); // 1MB page size
// void* a = alloc.New(15);  // Allocate 15 bytes
// void* b = alloc.New(27000);  // Allocate 27000 bytes
// void* c = alloc.New(150);  // Allocate 150 bytes
// void* d = alloc.New(150);  // Allocate 150 bytes
// void* e = alloc.New(150);  // Allocate 150 bytes
// alloc.Release(c); // Release a and b, but preserve c, d, e.
// alloc.Unalloc(e); // Free from e onward. (c, d still intact)

#ifndef FIFO_ALLOC_H
#define FIFO_ALLOC_H

//#define MEM_STATS 1

#include<assert.h>
#include<stdio.h>

// Used by FifoAllocator; probably not useful alone. See .cc file.
class FifoPage {
 public:
  explicit FifoPage(size_t size)
    : size_(size),
    page_((char*) GetPage(size)),
    end_(page_ + size_),
    next_(this),
    end_used_(page_),
    last_amt_(0) {
  }

  ~FifoPage() { DeletePage(page_, size_); }

  void Clear() { end_used_ = page_; }
  bool Empty() const { return end_used_ == page_; }

  void* New(size_t amount) {
    assert(amount >= 0);
    if (end_used_ + amount > end_)
      return NULL;
    void* pos = end_used_;
    end_used_ += amount;
    last_amt_ = amount;
    return pos;
  }

  void Unalloc(size_t amount) {
    assert(amount <= last_amt_);
    last_amt_ -= amount;
    end_used_ -= amount;
  }

  void Unalloc(void* pos) {
    assert(pos <= end_used_);
    Unalloc(end_used_ - (char*) pos);
  }

  void InsertPage(FifoPage* succ) {
    succ->next_ = next_;
    next_ = succ;
  }

  FifoPage* Next() { return next_; }

  bool InPage(void* ptr) {
    return ((char*) ptr >= page_ && (char*) ptr < end_used_);
  }

  bool AtEnd(void* ptr) {
    return ((char*) ptr == end_used_);
  }

  // For debugging
  void Show() {
    fprintf(stderr, "page[%p, )%p %p]", page_, end_used_, end_);
  }

 private:
  size_t size_;
  char* page_;
  char* end_;
  FifoPage* next_;
  char* end_used_;
  size_t last_amt_;

  static void* GetPage(size_t size);
  static void DeletePage(void* page, size_t size);
};


class FifoAllocator {
 public:
  explicit FifoAllocator(size_t page_size) : page_size_(page_size) {
    current_page_ = new FifoPage(page_size_);
    first_page_ = current_page_;
  }

  ~FifoAllocator();

  void* New(size_t amount) {
#ifdef MEM_STATS
    total_ += amount;
#endif
    void* result = current_page_->New(amount);
    if (result != NULL)
      return result;

    // No more room in current page, get another page.
    return FallbackNew(amount);
  }

  void Release(void* first_used);
  void ReleaseAll();

  void Unalloc(size_t amount) { current_page_->Unalloc(amount); }
  void Unalloc(void* pos) { current_page_->Unalloc(pos); }

  void Show();

#ifdef MEM_STATS
  static size_t Total() { return total_; }
#endif

 private:
  // Called when current page full
  void* FallbackNew(size_t amount);

  size_t page_size_;
  FifoPage* first_page_;
  FifoPage* current_page_;

#ifdef MEM_STATS
  static size_t total_;
#endif
};

#endif // FIFO_ALLOC_H
