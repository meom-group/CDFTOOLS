cmake_minimum_required(VERSION 3.2)
project(CDFTOOLS)
enable_language (Fortran)

set(EXECUTABLE_OUTPUT_PATH bin/)

string(REPLACE ":" ";" INCLUDE_DIRS $ENV{CPATH})
string(REPLACE ":" ";" LIB_DIRS $ENV{LD_LIBRARY_PATH})
include_directories(${INCLUDE_DIRS})
link_directories(${LIB_DIRS})

SET(GFORTRAN_FLAGS "-fno-second-underscore -ffree-line-length-256")

SET( CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${GFORTRAN_FLAGS}" )

# Modules

link_libraries(-lnetcdf -lnetcdff)

add_library(modcdfnames modcdfnames.f90)
add_library(modpoly modpoly.f90)

add_library(cdfio cdfio.F90)
target_link_libraries (cdfio modcdfnames )

add_library(cdftools cdftools.f90)
target_link_libraries (cdftools cdfio)

add_library(eos eos.f90)
target_link_libraries (eos cdfio)

add_library(modutils modutils.f90)
target_link_libraries (modutils cdfio)

# Executables

link_libraries(cdfio modcdfnames)

add_executable(cdfmoy cdfmoy.f90)
add_executable(cdfmoyt cdfmoyt.f90)
add_executable(cdfenstat cdfenstat.f90)
add_executable(cdfmoy_freq cdfmoy_freq.f90)
add_executable(cdfmoyuvwt cdfmoyuvwt.f90)
add_executable(cdfstd cdfstd.f90)
add_executable(cdfmoy_weighted cdfmoy_weighted.f90)
add_executable(cdfeke cdfeke.f90)
add_executable(cdfrmsssh cdfrmsssh.f90)
add_executable(cdfstdevw cdfstdevw.f90)
add_executable(cdfstdevts cdfstdevts.f90)
add_executable(cdfvT cdfvT.f90)
target_link_libraries (cdfvT modutils)
add_executable(cdfuv cdfuv.f90)
target_link_libraries (cdfuv modutils)
add_executable(cdfvsig cdfvsig.f90)
target_link_libraries (cdfvsig eos modutils)
add_executable(cdfspeed cdfspeed.f90)
add_executable(cdfimprovechk cdfimprovechk.f90)
add_executable(cdfstats cdfstats.f90)
target_link_libraries (cdfstats modutils)
add_executable(cdflinreg cdflinreg.f90)

link_libraries(eos)

add_executable(cdfbn2 cdfbn2.f90)
add_executable(cdfrichardson cdfrichardson.f90)
add_executable(cdfsig0 cdfsig0.f90)
add_executable(cdfsigntr cdfsigntr.f90)
add_executable(cdfspice cdfspice.f90)
add_executable(cdfsigi cdfsigi.f90)
add_executable(cdfsiginsitu cdfsiginsitu.f90)




