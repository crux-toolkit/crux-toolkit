#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include<deque>
#include "fifo_alloc.h"

using namespace std;

#define CHECK(x) { if (!(x)) { cerr << "CHECK FAILED: " #x << endl; abort(); } }

int Rand(int end) {
  // strictly less than end
  return random() % end;
}

typedef unsigned char uchar;

struct Entry {
  uchar* ptr;
  uchar start_val;
  uchar len;
};

typedef deque<Entry> Deq;

void CheckDeq(Deq* deq) {
  // printf("deq size: %d ", deq->size());
  for (Deq::iterator i = deq->begin(); i != deq->end(); ++i) {
    for (uchar j = 0; j < i->len; ++j) {
      if (i->ptr[j] != uchar(i->start_val + j)) {
	fprintf(stdout, "ERROR: i->ptr[j] = %d != i->start_val + j = %d\n",
		int(i->ptr[j]), int(i->start_val + j));
	abort();
      }
    }
  }
}

void Test2(int page_size, int max_to_alloc, int num_cycles) {
  printf("Test(page_size=%d, max_to_alloc=%d, num_cycles=%d)\n",
	 page_size, max_to_alloc, num_cycles);
  FifoAllocator fifo_alloc(page_size);
  Deq deq;
  uchar count;
  for (int i = 0; i < num_cycles; ++i) {
    // allocate some
    int num_to_alloc = Rand(max_to_alloc + 1);
    // printf("alloc %d ", num_to_alloc);
    for (int j = 0; j < num_to_alloc; ++j) {
      Entry e;
      e.len = Rand(10);
      e.ptr = (uchar*) fifo_alloc.New(e.len);
      e.start_val = count;
      for (int k = 0; k < e.len; ++k)
	e.ptr[k] = count++;
      deq.push_back(e);
    }
    CheckDeq(&deq);
    // deallocate some
    int num_to_dealloc = Rand(deq.size()+1);
    // printf("dealloc %d ", num_to_dealloc);
    if (num_to_dealloc >= deq.size()) {
      deq.clear();
      fifo_alloc.ReleaseAll();
    } else {
      uchar* first_used = deq[num_to_dealloc].ptr;
      fifo_alloc.Release(first_used);
      deq.erase(deq.begin(), deq.begin() + num_to_dealloc);
    }
    CheckDeq(&deq);
  }
}

void Test(int page_size, int num_to_alloc, int num_to_release, int num_cycles) {
  printf("Test(page_size=%d, num_to_alloc=%d, num_to_release=%d, num_cycles=%d)\n",
	 page_size, num_to_alloc, num_to_release, num_cycles);
  FifoAllocator fifo_alloc(page_size);
  int* arr[num_to_alloc];

  int count = 0;
  for (int i = 0; i < num_cycles; ++i) {
    for (int j = 0; j < num_to_alloc; ++j) {
      arr[j] = (int*) fifo_alloc.New(sizeof(arr[0]));
      CHECK(arr[j] != NULL);
      *(arr[j]) = count+j;
    }
    for (int j = 0; j < num_to_alloc; ++j)
      CHECK(*(arr[j]) == count+j);
    if (num_to_release < num_to_alloc) {
      fifo_alloc.Release(arr[num_to_release]);
    } else {
      fifo_alloc.ReleaseAll();
    }
    for (int j = num_to_release; j < num_to_alloc; ++j)
      CHECK(*(arr[j]) == count+j);
    count += num_to_alloc;
  }
}

int main(int argc, char* argv[]) {
  Test(21, 10, 3, 1);
  Test(21, 10, 3, 100);
  Test(20, 10, 3, 1);
  Test(20, 10, 3, 100);
  Test(21, 10, 10, 1);
  Test(21, 10, 10, 100);
  Test(20, 10, 10, 1);
  Test(20, 10, 10, 100);
  Test2(10000, 100, 1000);
  cout << "Passed" << endl;
  return 0;
}
