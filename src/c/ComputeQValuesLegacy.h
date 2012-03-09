/**
 *
 */
#ifndef ComputeQValuesLegacy_H
#define ComputeQValuesLegacy_H

#include "ComputeQValues.h"

class ComputeQValuesLegacy : public ComputeQValues {

 public:
  /**
   * \returns A blank ComputeQValues object.
   */
  ComputeQValuesLegacy(){};

  /**
   * Destructor
   */
  ~ComputeQValuesLegacy(){};

  /**
   * \returns The command name for ComputeQValues.
   */
  virtual std::string getName(){
    return "compute-q-values";
  };
  /**
   * Exclude this application from the usage statement.
   */
  virtual bool hidden(){
    return true;
  };
};



#endif //  ComputeQValuesLegacy_H


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
