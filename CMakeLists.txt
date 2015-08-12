cmake_minimum_required(VERSION 2.8.12)

PROJECT(vander)

SET(CMAKE_BUILD_TYPE RelWithDebInfo)

SET(MY_DIR ${vander_SOURCE_DIR})

find_package(PkgConfig)
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

add_definitions(${EIGEN_DEFINITIONS})

OPTION(ENABLE_CXX11 "Enable C++11 support" ON)

IF(ENABLE_CXX11)
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  IF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  ELSEIF(COMPILER_SUPPORTS_CXX0X)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  ELSE()
    MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  ENDIF()
ENDIF(ENABLE_CXX11)

INCLUDE_DIRECTORIES(EIGEN3_INCLUDE_DIRS
                AND "${MY_DIR}/include")


message(STATUS "        vander include path:\n")
message(STATUS               "${MY_DIR}/include")
message(STATUS "\neigen3 include path: ${EIGEN3_INCLUDE_DIRS}\n")

#Declare the executable target built from your sources
add_executable(vander vander.cpp)
