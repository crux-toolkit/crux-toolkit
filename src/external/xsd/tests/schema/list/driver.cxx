#include "test.hxx"

typedef xmlns::test::IntList<void, int> IntListImpl;
typedef xmlns::test::IntList<void, void> IntListVoidImpl;
typedef xmlns::test::IntComplex<void, int, char*> IntComplexImpl;

int
main ()
{
  IntListImpl int_list_impl;
  IntListVoidImpl int_list_void_impl;
  IntComplexImpl int_complex_impl;
}
