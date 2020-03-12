#!/usr/bin/env bash

my_name=`echo "$0" |sed 's,.*[\\/],,'`
if test $# -lt 1; then
    cat >&2 <<EOF

Usage: $my_name [nbThreads]

$my_name submit all jobs (one by one)

Options:
  --threads                       default: 1
  --batch                         default: ''
  --queue                         default: ''

EOF
    exit 1
fi

nbThreads=1
batch=''
memory=16000

while test $# -gt 0; do

    case "$1" in
   --threads) shift && nbThreads=`echo $1` ;;
     --batch) shift && batch=`echo $1` ;;
    --runMem) shift && memory=`echo $1` ;;
           *) # concatenate everything else to other_args
                # and hope that the user knows what they are doing.
               cat >&2 <<EOF
Error. Option $1 not supported.
Exiting...
EOF
               exit 1
               ;;

    esac
    shift

done


for parameterFile in $(find . -name "forwardMC*run?.txt")
do
  echo "Submitting job $parameterFile"	
  runTOPAS $parameterFile t=${nbThreads} batch=${batch} runMem=${memory}
done
