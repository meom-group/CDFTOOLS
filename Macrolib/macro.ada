# macro.ada for ada at IDRIS
#  $Rev$
#  $Date$
#  $Id$
# -------------------------------------------------------------
#
NCDF=/smplocal/pub/NetCDF/4.1.3
HDF5=/smplocal/pub/HDF5/1.8.9/seq
 
F90=ifort
MPF90=/opt/ibmhpc/pe1209/ppe.poe/bin/mpif90 -compiler ifort

# flag static is used to allow the use of CDFTOOLS in parallel with mpi_metamon
#FFLAGS= -O  $(NCDF) -assume byterecl -convert big_endian -CB -fpe0 -g -traceback -ftrapuv

FFLAGS= -O2 -I$(NCDF)/include -L$(NCDF)/lib -L$(HDF5)/lib -assume byterecl -convert big_endian -lnetcdff -lnetcdf  -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

#FFLAGS= -static -O  $(NCDF) -assume byterecl -convert big_endian
#FFLAGS= -static -O  $(NCDF) -assume byterecl -convert big_endian -g -traceback -fpe0 -ftrapuv -CB
LMPI=-lmpich

INSTALL=$(WORKDIR)/bin
INSTALL_MAN=$(WORKDIR)/man
