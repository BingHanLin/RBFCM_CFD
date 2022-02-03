# https://rix0r.nl/blog/2015/08/13/cmake-guide/
# https://stackoverflow.com/questions/16398937/cmake-and-finding-other-projects-and-their-dependencies
cmake_minimum_required(VERSION 3.12)

project(mainProject)


if(CMAKE_BUILD_TYPE MATCHES Release)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/bin/Release" CACHE PATH "The path to use for make install" FORCE)
else()
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/bin/Debug" CACHE PATH "The path to use for make install" FORCE)
endif()

##
## PROJECT NAME
##
project(MAIN_CMAKELIST)
message("CMAKE_SOURCE_DIR= ${CMAKE_SOURCE_DIR}")
message("CMAKE_BINARY_DIR= ${CMAKE_BINARY_DIR}")
##
## CONFIGURATION
##
set(CMAKE_CXX_STANDARD 17)

##
## TARGET
## create target and add include path
##
add_subdirectory(${PROJECT_SOURCE_DIR}/src/thirdParty)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/RFBCMCFD)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/app)

if(CATCH2TEST)
    # catch2 test
    add_subdirectory(${PROJECT_SOURCE_DIR}/src/test)
endif()
