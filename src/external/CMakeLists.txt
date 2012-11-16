include(ExternalProject)

# This CMakeLists.txt file is responsible for obtaining and building
# the required Proteowizard libraries. This includes dependencies like
# Boost, and Proteowizard proper.

# Set up quickbuild arguments
# quickbuild.sh is the script provided by proteowizard 
# to bootstrap bjam and to run bjam for building the external
# libraries and Proteowizard components..
set(protwiz_args ./quickbuild.sh)
set(protwiz_args ${protwiz_args} --without-binary-msdata)
set(protwiz_args ${protwiz_args} --layout=system)
set(protwiz_args ${protwiz_args} --prefix=${CMAKE_CURRENT_SOURCE_DIR})
set(protwiz_args ${protwiz_args} libraries)
set(protwiz_args ${protwiz_args} pwiz/data/common//pwiz_data_common)
set(protwiz_args ${protwiz_args} pwiz/data/identdata//pwiz_data_identdata)
set(protwiz_args ${protwiz_args} pwiz/data/identdata//pwiz_data_identdata_version)
set(protwiz_args ${protwiz_args} pwiz/data/msdata//pwiz_data_msdata)
set(protwiz_args ${protwiz_args} pwiz/data/msdata//pwiz_data_msdata_version)
set(protwiz_args ${protwiz_args} pwiz/data/proteome//pwiz_data_proteome)
set(protwiz_args ${protwiz_args} pwiz/utility/chemistry//pwiz_utility_chemistry)
set(protwiz_args ${protwiz_args} pwiz/utility/minimxml//pwiz_utility_minimxml)
set(protwiz_args ${protwiz_args} pwiz/utility/misc//pwiz_utility_misc)

# We want to avoid checking out the entire Proteowizard repository
# because of its massive size. We check out just the top level directory, 
# and them use 'svn update' to obtain the needed files and directories.
ExternalProject_Add(
  ProteoWizard
  SVN_REPOSITORY http://proteowizard.svn.sourceforge.net/svnroot/proteowizard/trunk/pwiz/
  SVN_REVISION -r 3950 --depth empty
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${protwiz_args}
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_links.sh --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)

# Since we're not checking out the entire Proteowizard tree,
# We need a separate step to use 'svn update' to pick up the files
# and directories we need.
ExternalProject_Add_Step(
  ProteoWizard
  update_selected_directories
  DEPENDEES download
  DEPENDERS configure
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/src/ProteoWizard
  COMMAND svn update --force -r 3950 libraries pwiz Jamroot.jam quickbuild.sh clean.sh
)
