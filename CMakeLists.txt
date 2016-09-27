cmake_minimum_required(VERSION 2.6)
project(CDFTOOLS)
enable_language (Fortran)


set(NETCDF_PATH "$ENV{EBROOTNETCDFMINFORTRAN}")

set(EXECUTABLE_OUTPUT_PATH bin/)

include_directories(${NETCDF_PATH}/include)
link_directories(${NETCDF_PATH}/lib64})

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

