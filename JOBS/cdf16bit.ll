#!/bin/ksh
# @ cpu_limit  = 7200
# @ data_limit = 1gb
# Nom du travail LoadLeveler
# @ job_name   = cdf16bit
# Fichier de sortie standard du travail
# @ output     = $(job_name).$(jobid)
# Fichier de sortie d'erreur du travail
# @ error      =  $(job_name).$(jobid)
# @ queue


CONFIG=ORCA025
CASE=G70
YEAR=1980

CONFCASE=${CONFIG}-${CASE}

for f in *gridU.nc ; do
  cdf16bit $f -check >> log
  mv cdf16bit.nc ${f}16
  mfput ${f}16 $CONFIG/$CONFCASE-MEAN16/$YEAR/
  mfput log $CONFIG/$CONFCASE-MEAN16/$YEAR/log.gridU
done
