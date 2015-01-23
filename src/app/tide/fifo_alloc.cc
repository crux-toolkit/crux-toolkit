// Benjamin Diament
//
// A FifoPage represents a large block of memory allocated at the system level.
// We currently use mmap and munmap to get system memory.
// FifoPage also contains a pointer to "next page" which induces a circular linked list
// (i.e. the last page points to the first).
// Schematically, we have the following,
//
//   FreePage1-->FreePage2-->...-->FreePageM-->UsedPage1-->UsedPage2-->...-->UsedPageN--+
//     ^                                                                                |
//     +-------------------------------------<------------------------------------------+
//
// current_page_ points to UsedPageN, and first_page_ points to UsedPage1.
//
// A system call to obtain new memory will occur whenever total allocations
// exceed a page's worth. 
// Pages become available for reuse when all contents are Release()'d.
//
// On Linux we use mmap to allocate memory and we mark the page as executable
// to provide run-time compilation of dot product calculations.

#include <sys/types.h>
#ifdef _MSC_VER
#include "mman.h"
#else
#include <sys/mman.h>
#endif
#include<stdlib.h>
#include<assert.h>
#include<iostream>
#include "fifo_alloc.h"

using namespace std;

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#ifdef MMAP_SENTINEL_CHECK
#undef NASSERT
#define SENTINEL_DATA_SIZE 100
#define SENTINEL_VALUE 0xAF

void FillSentinel(void* p, size_t size) {
  memset(p, SENTINEL_VALUE, size);
}

void CheckSentinel(void* p, size_t size) {
  for (size_t i = 0; i < size; ++i)
    CHECK(((char *) p)[i] == (char) SENTINEL_VALUE);
}

void* FifoPage::GetPage(size_t size) {
  // protections to allow exec (see above)
  int mmap_prot_mode = PROT_READ | PROT_WRITE | PROT_EXEC;
  // for sentinel data before and after
  size_t size_with_sentinels = size + 2 * SENTINEL_DATA_SIZE;
  void* p = mmap(0, size_with_sentinels, mmap_prot_mode, 
                 MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (p == NULL) {
    cerr << "Failed to allocate FifoPage of size " << size << ". Aborting\n";
    abort();
  }
  FillSentinel(p, SENTINEL_DATA_SIZE);
  FillSentinel((char *) p + size + SENTINEL_DATA_SIZE, SENTINEL_DATA_SIZE);
  void* tmp = (char *) p + SENTINEL_DATA_SIZE;
  // return (char *) p + SENTINEL_DATA_SIZE;
  cerr << "mmap'ed page at " << tmp << endl;
  return tmp;
}

void FifoPage::DeletePage(void* page, size_t size) {
  CheckSentinel((char *) page - SENTINEL_DATA_SIZE, SENTINEL_DATA_SIZE);
  CheckSentinel((char *) page + size, SENTINEL_DATA_SIZE);
  cerr << "munmap'ed " << page << endl;
  munmap((char *) page - SENTINEL_DATA_SIZE, size + 2 * SENTINEL_DATA_SIZE);
}
#else // MMAP_SENTINEL_CHECK
void* FifoPage::GetPage(size_t size) {
  // protections to allow exec (see above)
  int mmap_prot_mode = PROT_READ | PROT_WRITE | PROT_EXEC;
  void* p = mmap(0, size, mmap_prot_mode, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  if (p == NULL) {
    cerr << "Failed to allocate FifoPage of size " << size << ". Aborting\n";
    abort();
  }
  return p;
}

void FifoPage::DeletePage(void* page, size_t size) {
  munmap(page, size);
}
#endif // MMAP_SENTINEL_CHECK

void* FifoAllocator::FallbackNew(size_t amount) {    
  // Check if a free page is already in our linked list.
  FifoPage* free_page = current_page_->Next(); 
  if (free_page == first_page_) {  // No free page in linked list
    FifoPage* new_page = new FifoPage(page_size_);
    current_page_->InsertPage(new_page);
    current_page_ = new_page;
  } else {
    current_page_ = free_page;
  }
  assert(current_page_->Empty());

  if (amount > page_size_) { // CONSIDER: eliminate this restriction.
    cerr << "Requested " << amount << " bytes from FifoAllocator, "
	 << "but page size is " << page_size_ << ". Aborting\n";
    abort();
    return NULL;
  }

  void* result = current_page_->New(amount);
  assert(result != NULL);
  return result;
}

void FifoAllocator::Release(void* first_used) {
  // Release everything up to, but not including, first_used.
  while (!first_page_->InPage(first_used)) {
    if (first_page_ == current_page_) {
      if (first_page_->AtEnd(first_used)) {
	first_page_->Clear();
	return;
      }
      // Show();
      cerr << "Attempted release of address not in FifoAllocator range.\n";
      abort();
    }
    first_page_->Clear();
    first_page_ = first_page_->Next();
  }
  if (first_page_->AtEnd(first_used)) {
    first_page_->Clear();
    first_page_ = first_page_->Next();
  }
}

void FifoAllocator::ReleaseAll() {
  while (true) {
    first_page_->Clear();
    if (first_page_ == current_page_)
      break;
    first_page_ = first_page_->Next();
  }
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
  FifoPage* page = first_page_;
  while(true) {
    page->Show();
    if (page == current_page_)
      break;
    fprintf(stderr, " --> ");
    page = page->Next();
  }
  fprintf(stderr, "\n");
}

#ifdef MEM_STATS
size_t FifoAllocator::total_ = 0;
#endif
