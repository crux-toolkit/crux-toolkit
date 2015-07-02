#ifndef FIFO_ALLOC_H
#define FIFO_ALLOC_H

#include<assert.h>
#include<stdio.h>

class FifoPage {
 public:
  FifoPage(size_t size) : size_(size) {
    page_ = (char*) GetPage(size);
    end_ = page_ + size;
    next_ = this;
    ReleaseAll();
  }

  ~FifoPage() { DeletePage(page_, size_); }

  void InsertPage(FifoPage* succ) {
    succ->next_ = next_;
    next_ = succ;
  }

  FifoPage* Next() { return next_; }

  void* New(size_t amount);

  void Release(void* first_used) {
    assert(InRange(first_used));
    begin_used_ = (char*) first_used;
  }

  void ReleaseAll() {
    fprintf(stderr, "Release All on Page at %x\n", page_);
    end_used_ = page_;
    begin_used_ = NULL;
  }

  bool InPage(void* ptr) {
    return ((char*) ptr >= page_ && (char*) ptr < end_);
  }

  bool InRange(void* ptr) {
    if (!InPage((char*) ptr))
      return false;
    if (begin_used_ == NULL)
      return false;
    if (begin_used_ < end_used_)
      return ((char*) ptr >= begin_used_ && (char*) ptr < end_used_);
    return ((char*) ptr >= begin_used_ || (char*) ptr < end_used_);
  }

  void Show() {
    fprintf(stderr, "{ page(%x, %x) used(%x, %x) }",
	    page_, end_, begin_used_, end_used_);
  }

 private:
  size_t size_;
  char* page_;
  char* end_;
  FifoPage* next_;
  char* begin_used_; // beginning of used region; oldest allocated item. NULL if empty.
  char* end_used_; // just past used region; new allocs here.

  static void* GetPage(size_t size);
  static void DeletePage(void* page, size_t size);
};

class FifoAllocator {
 public:
  FifoAllocator(size_t page_size) : page_size_(page_size) {
    current_page_ = new FifoPage(page_size_);
  }

  ~FifoAllocator();

  void* New(size_t amount);
  void Release(void* first_used);
  void ReleaseAll();

  void Show();

 private:
  size_t page_size_;
  FifoPage* current_page_;
};

#endif // FIFO_ALLOC_H
