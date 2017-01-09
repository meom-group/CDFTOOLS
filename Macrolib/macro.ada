# macro.ada for ada at IDRIS
#  $Rev$
#  $Date$
#  $Id$
# -------------------------------------------------------------
#
 
F90=ifort
MPF90=mpiifort

# assume compiler wrapper is working for netcdf on ada 

# flag static is used to allow the use of CDFTOOLS in parallel with mpi_metamon
#FFLAGS= -O   -assume byterecl -convert big_endian -CB -fpe0 -g -traceback -ftrapuv

#FFLAGS= -O2 -assume byterecl -convert big_endian 
FFLAGS= -O2 -assume byterecl -convert big_endian -fpe0 -CB -ftrapuv -traceback -g -fp-model precise -fp-model source
#FFLAGS=-O2 -xSSE4.1 -ip -ftz -fpe3 -fno-alias -sox -assume byterecl -convert big_endian -fp-model precise -fp-model source

#OMP = -openmp
OMP = 

INSTALL=$(WORKDIR)/bin
INSTALL_MAN=$(WORKDIR)/man
