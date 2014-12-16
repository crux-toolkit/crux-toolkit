#include "HTMLWriter.h"

using namespace Crux;
using namespace boost;
using namespace std;


HTMLWriter::HTMLWriter() : PMCDelimitedFileWriter() {
  setWriteHTML(true);
}

HTMLWriter::~HTMLWriter() {
  closeFile();
}
