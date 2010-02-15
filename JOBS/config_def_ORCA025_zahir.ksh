#!/bin/ksh

# this is the config file for ORCA025 on zahir


#  $Rev$
#  $Date$
#  $Id$


# Name of CONFIG and CASE
CONFIG=ORCA025
CASE=G70

# check gaya .rhosts file !!!
# these variables are used only in the generic functions defined in Part II
#
LOCAL_SAVE=0  # set to 1 for using local storage instead of gaya for monitoring
USER=rcli300 ;   REMOTE_USER=rcli002  ; PREF=/u/rech/cli/$REMOTE_USER # PREF is the home of REMOTE_USER on remote machine

# Directory with the CDFTOOLS executable
CDFTOOLS=~rcli002/CDFTOOLS-2.1/

