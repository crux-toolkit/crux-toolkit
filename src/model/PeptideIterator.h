/**
 * \file PeptideIterator.h
 * Abstract class for returning peptides from a peptide source.
 */
#ifndef PETPIDE_ITERATOR_H
#define PETPIDE_ITERATOR_H

#include "Peptide.h"

class PeptideIterator{

 private:
  bool has_next_;
  
 protected:
  Crux::Peptide* next_peptide_;
  
  /**
   * \brief Prepare the iterator to return the next peptide.  
   * Derived classes must implement this method, pointing
   * next_peptide_ to the peptide that will be returned next and
   * returning true.  When there are no more peptides, it must set
   * has_next_ to false, set next_peptide_ to NULL and return false.
   * \returns true if there are more peptides to fetch or false if we are done
   */
  virtual bool queueNextPeptide() = 0;
  
  /**
   * Call this method in the constructor to set the first peptide.
   */
  void initialize(){
    has_next_ = queueNextPeptide();
  }

 public:
  virtual ~PeptideIterator(){};

  /**
   * Return the next peptide in the iterator or NULL if none are left.
   */
  Crux::Peptide* next(){
    Crux::Peptide* return_me = next_peptide_;
    has_next_ = queueNextPeptide();
    return return_me;
  }
  
  /**
   * Return true if there are more peptides to return or false if not.
   */
  bool hasNext(){
    return has_next_;
  }
  
};

#endif //PETPIDE_ITERATOR_H
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


