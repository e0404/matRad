#!/bin/bash

if test "x${MATRAD_PATH}" = 'x'; then
    cat >&2 <<EOF
matRad environment not defined
Load it first!
EOF

    exit 1
fi

my_name=`echo "$0" |sed 's,.*[\\/],,'`
if test $# -lt 2; then
    cat >&2 <<EOF

Usage: $my_name --matrad matRad_workspace_data --topas MCdata [--suffix suffix]

$my_name import forward MC results into resultGUI
EOF
    exit 1
fi

MATRAD_WORKSPACE=''
topas_workspace='MCdata.mat'
topas_metadata='MCparam.mat'
suffix='withMCresults'

while test $# -gt 0; do

    case "$1" in

      --matrad) shift && MATRAD_WORKSPACE=`echo $1` ;;
       --topas) shift && topas_workspace=`echo $1` ;;
    --metadata) shift && topas_metadata=`echo $1` ;;
      --suffix) shift && suffix=`echo $1` ;;
             *) # concatenate everything else to other_args
                # and hope that the user knows what they are doing.
                ;;
    esac
    shift
done

MATRAD_WORKSPACE="$(dirname $(readlink -e ${MATRAD_WORKSPACE}))/$(basename ${MATRAD_WORKSPACE})"

if [ ! -e $MATRAD_WORKSPACE ]
then
  echo "MatRad workspace file $MATRAD_WORKSPACE not found!"
  exit 1
fi

if [ ! -e $topas_metadata ]
then
  echo "MC metadata file $topas_metadata not found!"
  exit 1
fi

if [ ! -e $topas_workspace ]
then
  echo "MC workspace file $topas_workspace not found!"
  exit 1
fi

octave-cli --no-gui <<EOF 

disp('Reading results from matRad...');
load('$MATRAD_WORKSPACE')

disp(['resultGUI fieldnames: '])
disp(fieldnames(resultGUI))

disp('Reading metadata from TOPAS...');
load('$topas_metadata')

disp('Reading results from TOPAS...');
load('$topas_workspace')
fieldnames_resultMC = fieldnames(resultMC)
disp(fieldnames_resultMC)

for ix=1:length(fieldnames_resultMC)
  printf("Adding %s to resultGUI\n",fieldnames_resultMC{ix})
  resultGUI.(fieldnames_resultMC{ix}) = resultMC.(fieldnames_resultMC{ix});
end

clear resultMC, fieldnames_resultMC

if exist('OCTAVE_VERSION','builtin');
  % OCTAVE
  save('-v7',[ 'workspace_' '$suffix' '.mat']);

else
  % MATLAB
  save(fileName);
end

EOF
