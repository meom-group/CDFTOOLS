cmake_minimum_required(VERSION 3.2)
project(CDFTOOLS)
subdirs(src)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/)
set(ARCHIVE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib/)


string(REPLACE ":" ";" INCLUDE_DIRS $ENV{CPATH})
string(REPLACE ":" ";" LIB_DIRS $ENV{LD_LIBRARY_PATH})
include_directories(${INCLUDE_DIRS})
link_directories(${LIB_DIRS})

SET(GFORTRAN_FLAGS "-fno-second-underscore -ffree-line-length-256")

SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${GFORTRAN_FLAGS}" )

