# This is a CMake script for installing ProteoWizard for use in building crux.
#
# We have to copy the BOOST and PWIZ libraries out of their 
# build directory into our build lib directory
file(GLOB_RECURSE static_library_list ./*.lib)
file(COPY ${static_library_list} DESTINATION ./lib)
file(GLOB_RECURSE dll_list ./build/src/ProteoWizard/build-nt-x86/pwiz/data/vendor_readers/*.dll)
file(COPY ${dll_list} DESTINATION ./lib)

# Remove the compiler signature from boost library names
file(GLOB boost_library_list ./lib/libboost*.lib)
foreach(boost_library ${boost_library_list})
  string(REGEX REPLACE -vc14[0-9]- - renamed_boost_library ${boost_library})
  file(RENAME ${boost_library} ${renamed_boost_library})
endforeach()
