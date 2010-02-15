#!/bin/ksh
## Multistep job:
## step 1 : retrieve the restarts and forcing, run the model (186 procs)
## step 2 : build nc file and save them on gaya ( 16 procs)
#########################################################################
##     title of the run
# @ job_name = mmm

##     Output listing location
# @ output = $(job_name)-$(step_name).$(jobid)
# @ error  = $(output)

# @ step_name = p1
# @ job_type = serial
# @ cpu_limit = 36000
# @ data_limit = 1.3gb
# @ executable = ./CDFMOY
# @ queue

# @ step_name = p2
## @ dependency = (p1 == 0)
# @ job_type = serial
# @ cpu_limit = 36000
# @ data_limit = 1.3gb
# @ executable = ./CDFVT
# @ queue

# @ step_name = p3
# @ dependency = (p1 == 0 && p2 == 0 )
# @ job_type = serial
# @ cpu_limit = 7600
# @ data_limit = 0.3gb
# @ executable = ./MONITOR
# @ queue
