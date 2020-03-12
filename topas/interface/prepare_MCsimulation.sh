#!/bin/bash

# Script to prepare files for MC forward simulation of matRad plan
# Author: Lucas Norberto Burigo
# Contact: l.burigo@dkfz.de / burigolucas@gmail.com

set -e

my_name=`echo "$0" |sed 's,.*[\\/],,'`
if test $# -lt 1; then
    cat >&2 <<EOF

Usage: $my_name matRad_workspace [OPTION]
$my_name prepare files for forward MC simulation

Options:
  --beamSetup       "generic", "HITPS" [default: "generic"]
  --MCcode          "TOPAS", [default: "TOPAS"]
  --simDir          path, [default: "forwardMCdata"]
  --label           label, [default: "matrad_plan"]
  --fracHistories   fraction of particle histories
  --productionCut   production cut in mm
  --PBS             add flag to simulate PBS
  --custom          add flag to include B-Field or extra scorers

EOF
    exit 1
fi

#export MATRAD_WORKSPACE=`pwd`"/$1"
export MATRAD_WORKSPACE="$(dirname $(readlink -e $1))/$(basename $1)"

MC_SIM_DIR='.'
MC_TRANSPORT_CODE="TOPAS"
MC_BEAM_SETUP="generic" # either Generic or HITPS

shift 1

while test $# -gt 0; do

    case "$1" in

    --beamSetup) shift && MC_BEAM_SETUP=`echo $1` ;;
       --MCcode) shift && MC_TRANSPORT_CODE=`echo $1` ;;
       --simDir) shift && MC_SIM_DIR=`echo $1` ;;
        --label) shift && export MC_SIM_LABEL=`echo $1` ;;
--fracHistories) shift && export MC_FRAC_HISTORIES=`echo $1` ;;
--productionCut) shift && export MC_PRODUCTION_CUT_MM=`echo $1` ;;
       --custom) shift && export MC_CUSTOM='1' ;;
          --PBS) shift && export MC_PBS=`echo $1` ;;
              *) # anything else
                cat >&2 <<EOF
Error. Option $1 not supported.
Exiting...
EOF
                exit 1
                ;;

    esac
    shift

done

#MATRAD_PATH='/home/topas/data/matRad'
export MC_SIM_DIR
export MC_TRANSPORT_CODE
export MC_BEAM_SETUP

if [ ! -d $MC_SIM_DIR ]
then
    cat >&2 <<EOF

The directory $MC_SIM_DIR does no exist. Please, assign the correct path for the simulations using the option:
--simdir pathForMCsimulations

EOF
    exit 1
fi

echo "Preparing files for MC forward calculation for the workspace:"
echo "  $MATRAD_WORKSPACE"

octave-cli --no-gui --exec-path ${MATRAD_PATH} --path ${MATRAD_PATH} ${MATRAD_PATH}/matRad_exportMCinputFiles.m

cat <<EOF
Files for MC forward simulation created successfully.
EOF

exit 0
