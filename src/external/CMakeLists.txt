include(ExternalProject)

# Set up quickbuild arguments
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

ExternalProject_Add(
  ProteoWizard
  SVN_REPOSITORY http://proteowizard.svn.sourceforge.net/svnroot/proteowizard/trunk/pwiz/
  SVN_REVISION -r 3950
  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/build
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${protwiz_args}
  INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/proteowizard/create_links.sh --prefix=${CMAKE_CURRENT_SOURCE_DIR}
  BUILD_IN_SOURCE 1
)

