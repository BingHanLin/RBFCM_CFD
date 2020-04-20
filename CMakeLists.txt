# https://rix0r.nl/blog/2015/08/13/cmake-guide/
# https://stackoverflow.com/questions/16398937/cmake-and-finding-other-projects-and-their-dependencies
cmake_minimum_required(VERSION 3.1)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

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
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/thirdParty)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/RFBCM_CFD)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/app)