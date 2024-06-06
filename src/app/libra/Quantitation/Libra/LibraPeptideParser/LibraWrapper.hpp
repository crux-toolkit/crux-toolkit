#ifndef LIBRAWRAPPER_HPP
#define LIBRAWRAPPER_HPP

#include "LibraConditionHandler.hpp"
#include "LibraResult.hpp"
#include "LibraSummary.hpp"
#include "Quantitation.hpp"

#include<vector>
#include<utility>
#include "typedefs.hpp"

class LibraWrapper
{
private:
  mzParser::RAMPFILE* pFI;
  
  char* mzXMLFile;
  
  LibraConditionHandler* pLibraConditionHandler;
  
public:
  
  LibraWrapper( LibraConditionHandler* , char*, mzParser::RAMPFILE* );
  
  LibraWrapper( char* , char*, mzParser::RAMPFILE* );
  
  ~LibraWrapper();
  
  LibraSummary* getLibraSummary();
  
  LibraResult* getLibraResult(int scanNum);
  
};

#endif
