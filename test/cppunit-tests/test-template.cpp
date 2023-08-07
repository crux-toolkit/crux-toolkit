// HOW TO USE THIS TEMPLATE:
// 1. Copy the file, renaming it TestMyClass.cpp where 'MyClass' is replaced
//    by the name of the new class being tested.
// 2. Add the filename to Makefile as part of the TEST variable.
// 3. Add the new file to SVN.
// 4. In the new file, delete this comment.
// 5. Replace all instances of [class] with the name of the new class.
// 6. Implement setUp() and tearDown().
// 7. Replace myTestName with the name of your first test.
// 8. Implement that test and others.
// SEE ALSO test-template.h

#include <cppunit/config/SourcePrefix.h>
#include "Test[class].h"
#include "parameter.h" 

CPPUNIT_TEST_SUITE_REGISTRATION( Test[class] );

void Test[class]::setUp(){
  // initialize_parameters();  // accessing any parameter values requires this

  // initialize variables to use for testing
}

void Test[class]::tearDown(){
  // delete anything you allocated
}

void Test[class]::myTestName(){
  // available tests include 
  CPPUNIT_ASSERT(b == a);

}










