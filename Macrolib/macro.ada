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

FFLAGS= -O2 -assume byterecl -convert big_endian 

INSTALL=$(WORKDIR)/bin
INSTALL_MAN=$(WORKDIR)/man
