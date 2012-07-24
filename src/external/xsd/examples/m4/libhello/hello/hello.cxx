#include <hello/hello.hxx>

#include <iostream>

void hello::
say (char const* phrase)
{
  std::cerr << phrase << std::endl;
}
