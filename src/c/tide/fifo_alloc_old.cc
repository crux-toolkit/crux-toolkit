#include <sys/types.h>
#include <sys/mman.h>
#include<stdlib.h>
#include<assert.h>
#include<iostream>
#include "fifo_alloc.h"

using namespace std;

void* FifoPage::New(size_t amount) {
  assert(amount >= 0);
  if (begin_used_ == NULL) {
    if (amount > size_)
      return NULL;
    begin_used_ = page_;
    end_used_ = page_ + amount;
    return page_;
  }
  if (begin_used_ < end_used_) {
    if (amount <= (end_ - end_used_)) {
      void* pos = end_used_;
      end_used_ += amount;
      return pos;
    }
    if (amount <= (begin_used_ - page_)) {
      end_used_ = page_ + amount;
      return page_;
    }
    return NULL;
  }
  if (amount > (begin_used_ - end_used_))
    return NULL;
  void* pos = end_used_;
  end_used_ += amount;
  return pos;
}

void* FifoPage::GetPage(size_t size) {
  // cerr << "Getting new page of size " << size << endl;
  void* p = mmap(0, size, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (p == NULL) {
    cerr << "Failed to allocate FifoPage of size " << size << ". Aborting\n";
    abort();
  }
  return p;
}

void FifoPage::DeletePage(void* page, size_t size) {
  munmap(page, size);
}

void* FifoAllocator::New(size_t amount) {
  void* result = current_page_->New(amount);
  if (result != NULL)
    return result;

  FifoPage* next = current_page_->Next();
  if (next != current_page_) {
    result = next->New(amount);
    if (result != NULL) {
      current_page_ = next;
      return result;
    }
  }

  if (amount > page_size_) {
    cerr << "Requested " << amount << " bytes from FifoAllocator, "
	 << "but page size is " << page_size_ << ". Aborting\n";
    abort();
    return NULL;
  }

  FifoPage* new_page = new FifoPage(page_size_);
  current_page_->InsertPage(new_page);
  current_page_ = new_page;

  result = current_page_->New(amount);
  assert(result != NULL);
  return result;
}

void FifoAllocator::Release(void* first_used) {
  FifoPage* page = current_page_->Next();
  while (!page->InPage(first_used)) {
    if (page == current_page_) {
      cerr << "Attempted release of page not in FifoAllocator range.\n";
      abort();
    }
    page->ReleaseAll();
    page = page->Next();
  }
  if (!page->InRange(first_used)) {
    cerr << "Attempted release of page not in FifoAllocator range.\n";
    abort();
  }
  page->Release(first_used);
}

void FifoAllocator::ReleaseAll() {
  FifoPage* page = current_page_;
  do {
    page = page->Next();
    page->ReleaseAll();
  } while (page != current_page_);
}

FifoAllocator::~FifoAllocator() {
  FifoPage* page = current_page_;
  do {
    FifoPage* next = page->Next();
    delete page;
    page = next;
  } while (page != current_page_);
}

void FifoAllocator::Show() {
  FifoPage* page = current_page_;
  while(true) {
    page = page->Next();
    page->Show();
    if (page == current_page_)
      break;
    fprintf(stderr, " --> ");
  }
  fprintf(stderr, "\n");
}
