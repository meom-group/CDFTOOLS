# macro.ada for ada at IDRIS
#  $Rev$
#  $Date$
#  $Id$
# -------------------------------------------------------------
#
 
F90=ifort
MPF90=mpiifort

# assume compiler wrapper is working for netcdf on ada 

#OMP = -openmp
OMP = 
#NC4 = -D key_netcdf4
NC4 = 

#CMIP6 = -D key_CMIP6
CMIP6 = 

# flag static is used to allow the use of CDFTOOLS in parallel with mpi_metamon
#FFLAGS= -O  $(NC4) $(CMIP6)  -assume byterecl -convert big_endian -CB -fpe0 -g -traceback -ftrapuv $(OMP)

#FFLAGS= -O2 $(NC4) $(CMIP6) -assume byterecl -convert big_endian  $(OMP)
FFLAGS= -O2 $(NC4) $(CMIP6) -assume byterecl -convert big_endian -fpe0 -CB -ftrapuv -traceback -g -fp-model precise -fp-model source $(OMP)
#FFLAGS=-O2 $(NC4) $(CMIP6) -xSSE4.1 -ip -ftz -fpe3 -fno-alias -sox -assume byterecl -convert big_endian -fp-model precise -fp-model source $(OMP)

INSTALL=$(WORKDIR)/bin
INSTALL_MAN=$(WORKDIR)/man
