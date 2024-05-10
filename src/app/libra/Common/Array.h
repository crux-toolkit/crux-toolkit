#ifndef ARRAY_INCLUDED 
#define ARRAY_INCLUDED

// copied shamelessly from Lippman
// chrisa 18-Nov-94

//
// MODIFIED for speed, uses bitwise copy instead of assignment 8-22-05 BSP
//
// if you need a container class that employs constructors and destructors, use STL::vector
// but this is great for pointers, ints, doubles, structs etc - but don't expect constructors
// or destructors to fire!!!!
//
// There are a handful of wrapper classes at the bottom of this file for handling arrays
// of things that do need deletion.
//


#ifdef SCCS_ID
static const char* arraySccsId = "@(#)array.h	1.4 09/15/95 11:43:21";
#endif

#include <iostream>

#ifdef _DEBUG
//#define ARRAY_H_DEBUG
#endif

#ifdef ARRAY_H_DEBUG 
#include "assert.h"
#define ARRAY_H_DEBUG_ASSERT assert
#define ARRAY_H_DEBUG_STATEMENT(x) x 
#else
#define ARRAY_H_DEBUG_ASSERT(x) // nothing
#define ARRAY_H_DEBUG_STATEMENT(x) // nothing 
#endif

#include "constants.h" // defines strCopy

const int iDefArrayCapacity[2] = {5, 100}; // we'll grow in order 0,5,100,200,400,800,...

template<class T>
class Array {
public:
  // ctor with no args gets empty array.  default ctor
  // for Arrays that are components of other classes.
  Array();

  // ctor with one int arg allocates but does not initialize
  // elements - so size() returns iS but this->[iS-1] is uninit.
  Array(const int iS);

  // ctor can initialize from array of T and int Capacity
  Array(const T *tyArr, const int iCapacity_);

  // preset capacity to avoid grow calls - leaves iUsed__ at 0
  void reserve(const int capacity);

  // set array to a particular size.  deletes prev contents!
  // useful if you want to go on and populate element[n] directly,
  // but note that insertAtEnd will create an iS+1th element
  void clearAndSetSize(const int iS);

  // set array to a particular size.  does not delete prev contents
  // if iS is greater than current size, end elements are uninit
  void setSize(const int iS);
    
  // deletes prev contents!
  void clear() {
     clearAndSetSize(0);
  }
    
  // copy constructor.  size may change, since only iUsed__ are copied.
  Array(const Array<T> &arrRef);

  ~Array();

  Array<T>& operator=(const Array<T>& aRef);

  // returns the number of used elements - not the same as the storage capacity
  // remove any excess capacity
  void trim();
  inline int size() const {
     return iUsed_; 
  }
  inline int length() const { // synonym for size() (it's evil to have synonyms! but changing code is bad too)
     return iUsed_; 
  }
  
  // wipe memory
  void nullify() {
     memset(arrPtr_,0,iCapacity_*sizeof(T));
  }
  
  // insert a new element of T at pos
  // will bomb if you attempt to insert past first free
  // (i.e. iPos > size() )
  void insert(const int iPos, const T& newElement);

  // insert at the end of the array
  void insertAtEnd(const T& newElement);

  // move nth element to end
  void moveLast(int iPos) {
     if ((iPos>=0)&&(iPos<iUsed_-1)) {
        T tmp = this->arrPtr_[iPos];
        memmove(this->arrPtr_+iPos,this->arrPtr_+iPos+1,sizeof(T)*(iUsed_-(iPos+1)));
        this->arrPtr_[iUsed_-1] = tmp;
     }
  }

  inline void replace(const int iPos, const T& newElement) {
    arrPtr_[iPos] = newElement;
  }

  // const element access
  inline const T& operator [] (int iIndex) const { 
     ARRAY_H_DEBUG_ASSERT((iIndex >= 0) && (iIndex < iUsed_));
     return this->arrPtr_[iIndex];
  }
  
  // non-const element access
  inline T& operator [] (int iIndex) { 
     ARRAY_H_DEBUG_ASSERT((iIndex >= 0) && (iIndex < iUsed_));
     return this->arrPtr_[iIndex];
  }

  void remove(const int iPos);

  void removeByValue(const T& killMe);

  int findByValue(const T& findMe) const; 

  // reverses the order of elements in the array, i.e.
  // a[0] <swapped with> a[iUsed_], etc.
  void reverseOrder();

  // sort according to user supplied compare function
  void sort(int (*compare )(const void *elem1, const void *elem2 ) ) {
     qsort( arrPtr_, iUsed_, sizeof(T), compare );
  }

  // Binary search, only useful on sorted data
  T* find_sorted(const T* what,int (*compare )(const void *elem1, const void *elem2 ) ) {
     return arrPtr_?(T*) bsearch(what, arrPtr_, iUsed_, sizeof(T), compare ):NULL; 
  }

protected:
  // init the array of T with passed array
  void init (const T* array, const int iU);

  // allocs new mem of (factor * current) Capacity, 
  // copies old to new.
  void grow();
  T *arrPtr_;     // ptr to array of T 
  int iCapacity_;        // size of currently allocated array
  int iUsed_;        // how much is used (iUsed-1 == <index last elem>)
};

template<class T>Array<T>::Array() : arrPtr_(0), iCapacity_(0), iUsed_(0) {
}

template<class T>Array<T>::Array(const int iS) : arrPtr_(0), iCapacity_(0), iUsed_(0) {
   clearAndSetSize(iS);  // though not initialized, they are available
}

template<class T>Array<T>::Array(const T *tyArr, const int iCapacity_) : arrPtr_(0), iCapacity_(0), iUsed_(0) { 
  init(tyArr, iCapacity_); 
}

template<class T> void Array<T>::clearAndSetSize(const int iS) {
   if (iS > iCapacity_) {// prev contents lost
      free(arrPtr_); 
      arrPtr_ = (T*)malloc(sizeof(T)*(iCapacity_ = iS));
   }
    iUsed_ = iS;  // allow refs to these (empty) elems
}

template<class T> void Array<T>::setSize(const int iS) {
   arrPtr_ = (T*)realloc(arrPtr_,sizeof(T)*(iCapacity_ = iS));
   iUsed_ = iS;  // allow refs to these (possibly empty) elems
}

template<class T> void Array<T>::reserve(const int iS) { // preset capacity to avoid grow
    // prev contents lost
    clearAndSetSize(iS);
    iUsed_ = 0;  // nothing occupied yet
}

template<class T> Array<T>::Array(const Array<T> &arrRef) : arrPtr_(0), iCapacity_(0), iUsed_(0) { 
  init(arrRef.arrPtr_, arrRef.iUsed_); 
}

template<class T> Array<T>::~Array() { 
  free (this->arrPtr_);
}

template<class T> Array<T>& Array<T>::operator=(const Array<T>& aRef) {
    // never copy to yourself
    if (this != &aRef) {
      free(arrPtr_);      // get rid of the old storage
      arrPtr_ = NULL;
      // re-allocate and copy the T refs stuff into yours
      // note this ignores actual storage size of passed ref
      init(aRef.arrPtr_, aRef.iUsed_);
    }
    return *this;      // and return ref to yourself
}

template<class T> void Array<T>::insert(const int iPos, const T& newElement) {
    ARRAY_H_DEBUG_ASSERT(iPos <= iUsed_);
    ARRAY_H_DEBUG_STATEMENT(T hold = iUsed_?arrPtr_[iUsed_-1]:newElement);

    if ((iUsed_+1) > iCapacity_) {
       grow();  // no room? grow by default factor
    }
    // shift contents of array from iPos to end up 1
    // for (int i = iUsed_-1; i >= iPos; i--) arrPtr_[i+1] = arrPtr_[i];
    memmove(arrPtr_+iPos+1,arrPtr_+iPos,(iUsed_-iPos)*sizeof(T));
    // assign elem at iPos to passed ref
    arrPtr_[iPos] = newElement;
    iUsed_++;  // inc current size count
    ARRAY_H_DEBUG_ASSERT(!memcmp(&hold,&arrPtr_[iUsed_-1],sizeof(T)));
}

template<class T> void Array<T>::insertAtEnd(const T& newElement) {
   if ((iUsed_+1) > iCapacity_) {
      grow();  // no room? grow by default factor
   }
    arrPtr_[iUsed_++] = newElement;
}

template<class T> void Array<T>::remove(const int iPos) {
  ARRAY_H_DEBUG_ASSERT((iPos >= 0) && (iPos < iUsed_));
  //for (int i = iPos; i < (iUsed_-1); i++)    arrPtr_[i] = arrPtr_[i+1];
  memmove(arrPtr_+iPos,arrPtr_+iPos+1,((iUsed_-iPos)-1)*sizeof(T));
  iUsed_--;  // dec current size count, storage unchanged
  ARRAY_H_DEBUG_ASSERT((iPos==iUsed_)||(!memcmp(&arrPtr_[iUsed_],&arrPtr_[iUsed_-1],sizeof(T))));
}

template<class T> void Array<T>::removeByValue(const T& killMe) {
    // read backward through the array
    for (int i = (iUsed_ - 1); i >= 0; i--) {
      if (arrPtr_[i] == killMe) {
        remove(i);
      }
    }
}

template<class T> int Array<T>::findByValue(const T& findMe) const {
    for (int i = 0; i < iUsed_; i++) {
      if (arrPtr_[i] == findMe) { 
         return(i); 
    }
    }
    return (-1);
}


template<class T> void Array<T>::reverseOrder() {
    if (iUsed_ < 2) return;
    for (int i = 0; i < (iUsed_ / 2); i++) {
      T tTemp(arrPtr_[i]);
      arrPtr_[i] = arrPtr_[iUsed_ - i - 1];
      arrPtr_[iUsed_ -i -1] = tTemp;
    }
}

template<class T> void Array<T>::init (const T* array, const int iU) {
    // storage size is recalculated
    arrPtr_ = (T*)malloc(sizeof(T)*(iCapacity_ = iU));
    iUsed_ = iU;  // let them see what they passed
    ARRAY_H_DEBUG_ASSERT( arrPtr_ != 0 );
    //for (int i = 0; i < iCapacity_; i++) arrPtr_[i] = array[i];
    memmove(arrPtr_,array,iUsed_*sizeof(T));
}

template<class T> void Array<T>::grow() {
  // increase the storage capacity
  if (iCapacity_ == 0) {
    iCapacity_ = iDefArrayCapacity[0]; // start small
  } else if (iCapacity_ == iDefArrayCapacity[0]) {
    iCapacity_ = iDefArrayCapacity[1]; // jump from small to biggish
  } else { 
    iCapacity_ *= 2; // double up
  }
  arrPtr_ = (T*)realloc(arrPtr_,sizeof(T)*iCapacity_);  
  ARRAY_H_DEBUG_ASSERT(arrPtr_ != NULL);
}

// trim any excess capacity
template<class T> void Array<T>::trim() {
   if (iUsed_ < iCapacity_) {
      arrPtr_ = (T*)realloc(arrPtr_,sizeof(T)*(iCapacity_ = iUsed_));  ARRAY_H_DEBUG_ASSERT(arrPtr_ != 0);
   }
}

//
// convenience class for arrays of arrays
//
template<class T>
class ArrayArray : public Array<Array<T>*> {
public:
  ~ArrayArray() {
     for (int i=this->size();i--;) {
       free(this->arrPtr_[i]);
     }
   }
};

//
// convenience class for arrays of arrays of pointers that need freeing
//
template<class T>
class PointerArrayArray : public Array<Array<T>*> {
public:
   ~PointerArrayArray() {
      for (int i=this->size();i--;) { // gcc is only happy when we use this-> (why?)
         for (int j=this->arrPtr_[i]->size();j--;) {
	   free((*(this->arrPtr_[i]))[j]);
            // beware duplicate refs
            for (int ii=i;ii>=0;ii--) {
               for (int jj=this->arrPtr_[ii]->size();jj--;) {
                  if (((i!=ii)||(j!=jj)) &&
                     ((*(this->arrPtr_[i]))[j] == 
                      (*(this->arrPtr_[ii]))[jj])) {
                     (*(this->arrPtr_[ii]))[jj] = NULL;
                  }
               }
            }
         }
         free( this->arrPtr_[i] );
      }
   }
};

//
// convenience class for arrays of strings that need freeing
//
class StringArray : public Array<char *> {
public:
   ~StringArray() {
     for (int j=0;j<this->size();j++) {
       delete [] this->arrPtr_[j];
     }
   }
   StringArray &operator=(const StringArray& aRef) {
    // never copy to yourself
    if (this != &aRef) {
      for (int i=this->size();i--;) {
	free(this->arrPtr_[i]);
       this->arrPtr_[i] = NULL;
      }
      // re-allocate and copy the T refs stuff into yours
      // note this ignores actual storage size of passed ref
      reserve(aRef.iUsed_);
      for (int j=0;j<aRef.size();j++) {
         insertAtEnd(strCopy(aRef.arrPtr_[j]));
      }
    }
    return *this;      // and return ref to yourself
   }
   int findByStringValue(const char *findMe) const {
    for (int i = 0; i < iUsed_; i++) {
      if (strcmp(arrPtr_[i], findMe)==0) { 
         return(i); 
      }
    }
    return (-1);
}

};


//
// convenience class for arrays of arrays of strings that need freeing
//
class StringArrayArray : public Array<StringArray *> {
public:
   ~StringArrayArray() {
      for (int i=this->size();i--;) {
	free(this->arrPtr_[i]);
      }
   }
};

#endif
   
