# https://rix0r.nl/blog/2015/08/13/cmake-guide/
# https://stackoverflow.com/questions/16398937/cmake-and-finding-other-projects-and-their-dependencies
cmake_minimum_required(VERSION 3.1)

##
## PROJECT NAME
##
project(MAIN_CMAKELIST)

##
## CONFIGURATION
##

##
## TARGET
## create target and add include path
##
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/thirdParty)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/RFBCM_CFD)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/app)