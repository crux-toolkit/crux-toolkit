// Benjamin Diament

// Class to represent a sized array, but with constant capacity, for
// speed. No bounds checking or other checking is performed.
//
// In FixedCapacityArray size_ always represents the client's interpretation
// of the number of elements in the array while capacity supplied to Init()
// never changes.
// The vector<> template has the problem that if you shrink the size then the
// the excess capacity will be freed. If later the array grows an allocation
// is performed. I couldn't find another way around this inefficiency.
// 
// iterator supplied to match vector<> template usage.
//
// Init() will permit optional use of a FifoAllocator for allocation.

#ifndef FIXED_CAP_ARRAY_H
#define FIXED_CAP_ARRAY_H

#include "fifo_alloc.h"

template <class C>
class FixedCapacityArray {
 public:
  explicit FixedCapacityArray(int capacity)
    : data_(new C[capacity]), size_(0), del_(true) {
  }

  // must call Init before use
  FixedCapacityArray() 
    : data_(NULL), 
    size_(0),
    del_(true) {
  }

  void Init(int capacity) { data_ = new C[capacity]; }

  void Init(FifoAllocator* fifo_alloc, int capacity) {
    if (fifo_alloc == NULL) {
      Init(capacity);
      return;
    }
    void* buffer = fifo_alloc->New(capacity * sizeof(C));
    data_ = (C*) buffer;
    del_ = false;
  }

  ~FixedCapacityArray() { if (del_) delete[] data_; }
  
  void clear() { size_ = 0; }

  bool empty() const { return size_ == 0; }

  C operator[](int index) const { return data_[index]; }
  C& operator[](int index) { return data_[index]; }

  C* data() { return data_; }
  int size() const { return size_; }
  void set_size(int size) { size_ = size; }

  void push_back(const C& elt) { data_[size_++] = elt; }
  C back() { return data_[size_-1]; }

  typedef C* iterator;
  typedef const C* const_iterator;

  iterator begin() { return data_; }
  iterator end() { return data_ + size_; }
  const_iterator begin() const { return data_; }
  const_iterator end() const { return data_ + size_; }

 private:
  C* data_;
  int size_;

  bool del_; // True if new/delete used. False if FifoAllocator used, 
             // in which case client deallocates.  
};

#endif // FIXED_CAP_ARRAY_H
