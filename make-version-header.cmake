set(CRUX_VERSION_FULL "${CRUX_VERSION}")

if (EXISTS "${SOURCE_DIR}/.git/HEAD")
  file(READ "${SOURCE_DIR}/.git/HEAD" PROJECT_SOURCE_VERSION)
  if ("${PROJECT_SOURCE_VERSION}" MATCHES "^ref: ")
    string(REGEX REPLACE "^ref: (.*)[\r\n]" "\\1" PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
    if (EXISTS "${SOURCE_DIR}/.git/${PROJECT_GIT_REF}")
      file(READ "${SOURCE_DIR}/.git/${PROJECT_GIT_REF}" PROJECT_SOURCE_VERSION)
      string(STRIP "${PROJECT_SOURCE_VERSION}" PROJECT_SOURCE_VERSION)
      set(CRUX_VERSION_FULL "${CRUX_VERSION_FULL}-${PROJECT_SOURCE_VERSION}")
    else (EXISTS "${SOURCE_DIR}/.git/${PROJECT_GIT_REF}")
      MESSAGE("Version header error: ${SOURCE_DIR}/.git/${PROJECT_GIT_REF} not found")
    endif (EXISTS "${SOURCE_DIR}/.git/${PROJECT_GIT_REF}")
  else ("${PROJECT_SOURCE_VERSION}" MATCHES "^ref: ")
    MESSAGE("Version header error: ${SOURCE_DIR}/.git/HEAD does not match pattern")
  endif ("${PROJECT_SOURCE_VERSION}" MATCHES "^ref: ")
else (EXISTS "${SOURCE_DIR}/.git/HEAD")
  MESSAGE("Version header error: ${SOURCE_DIR}/.git/HEAD not found")
endif (EXISTS "${SOURCE_DIR}/.git/HEAD")

file(WRITE crux_version.h.txt "#define CRUX_VERSION \"${CRUX_VERSION_FULL}\"\n")

# Copy the file to the final header if the version changed
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different crux_version.h.txt crux_version.h)

