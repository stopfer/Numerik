cmake_minimum_required(VERSION 2.6)
 
project(Numerik)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")

if (CMAKE_BUILD_TYPE STREQUAL "")
  # CMake defaults to leaving CMAKE_BUILD_TYPE empty. This screws up
  # differentiation between debug and release builds.
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

#hdnum include dir
set(HDNUM_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/hdnum/")

#libgmp installed?
find_package(GMP)

#gnuplot installed?
find_package(Gnuplot REQUIRED)

# LaTeX installed?
FIND_PACKAGE(LATEX)

if(PDFLATEX_COMPILER)
    include(UseLATEX)
    set(LATEX_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/doc")
else(PDFLATEX_COMPILER)
    message("LaTeX not found!")
endif(PDFLATEX_COMPILER)

#add include directories
include_directories(${GMP_INCLUDE_DIR} ${HDNUM_INCLUDE_DIR})

add_subdirectory(src)