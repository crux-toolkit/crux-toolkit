// HOW TO USE THIS TEMPLATE:
// 1. Copy the file, renaming it TestMyClass.h where 'MyClass' is replaced
//    by the class name being tested.
// 2. In the new file, delete this comment
// 3. Replace all instances of [class] with the name of the new class.
// 4. Replace myTestName with the name of the first test in the .cpp file
// 5. Declare each new test.
// 6. Add additional CPPUNIT_TEST( name ); lines for each new test
// 7. Add the new file to SVN.
// SEE ALSO test-template.cpp

#ifndef CPP_UNIT_[class]_H
#define CPP_UNIT_[class]_H

#include <cppunit/extensions/HelperMacros.h>
#include "[class].h"

class Test[class] : public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE( Test[class] );
  CPPUNIT_TEST( myTestName );
  CPPUNIT_TEST_SUITE_END();
  
 protected:
  // variables to use in testing

 public:
  void setUp();
  void tearDown();

 protected:
  void myTestName();
};

#endif //CPP_UNIT_[class]_H
