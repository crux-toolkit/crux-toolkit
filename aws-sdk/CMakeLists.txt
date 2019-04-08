
include(ExternalProject)


set(BUILD_SHARED_LIBS 0)
cmake_minimum_required(VERSION 2.8.4)
cmake_policy(VERSION 2.8.4)

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/build/aws-deps")


foreach(DEPENDENCY_NAME aws-c-common aws-c-event-stream aws-checksums)
        ExternalProject_add(
            ${DEPENDENCY_NAME}
            PREFIX aws-sdk/build/${DEPENDENCY_NAME}
            INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
            GIT_REPOSITORY "https://github.com/awslabs/${DEPENDENCY_NAME}.git"
            GIT_PROGRESS 1
            CMAKE_ARGS
                "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
            UPDATE_COMMAND ""
        )
endforeach(DEPENDENCY_NAME)

set(LIB_NAME aws-sdk-cpp)
ExternalProject_add(
        ${LIB_NAME}
        PREFIX aws-sdk/build/${LIB_NAME}
        GIT_REPOSITORY "https://github.com/aws/${LIB_NAME}.git"
        GIT_PROGRESS 1
        CMAKE_ARGS
                "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                "-DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
                "-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}"
                "-DBUILD_ONLY=s3"
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
        UPDATE_COMMAND ""
        #       BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
        #       COMMAND "${CMAKE_MAKE_PROGRAM} install"
)

