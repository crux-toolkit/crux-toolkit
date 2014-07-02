if (EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  include(FindSubversion)

# extract working copy information for SOURCE_DIR into MY_XXX variables
  Subversion_WC_INFO(${SOURCE_DIR} MY)

  file(WRITE crux_version.h.txt "#define CRUX_VERSION \"${CRUX_VERSION}.${MY_WC_REVISION}\"\n")

# Copy the file to the final header only if the version changes
# reduces needless rebuilds
  execute_process(
    COMMAND 
      ${CMAKE_COMMAND} -E copy_if_different crux_version.h.txt crux_version.h
  )
else (EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  file(WRITE crux_version.h.txt "#define CRUX_VERSION \"${CRUX_VERSION}\"\n")

# Copy the file to the final header only if the version changes
# reduces needless rebuilds
  execute_process(
    COMMAND 
      ${CMAKE_COMMAND} -E copy_if_different crux_version.h.txt crux_version.h
  )
endif (EXISTS ${CMAKE_SOURCE_DIR}/.svn)
