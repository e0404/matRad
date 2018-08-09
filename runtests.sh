#!/usr/bin/env bash

echo 'what up'
## Make sure some failures are detected by the CI runners
function exitIfError {
	# pass "$?" as argument: i.e. the exit status of the last call
	if [ "$1" -ne 0 ]; then
		exit $1;
	fi
}

## Handle Runner and Switches variables
Runner=$1
Switches=$2
if [ -z "$Runner" ] ; then
	Runner="octave"
fi
if [ -z "$Switches" ] ; then
	case "$Runner" in
		*matlab* )
			Switches="-nodesktop -r"
			;;

		*octave* )
			Switches="--no-gui --eval"
			;;

		* )
			# Fall back to Octave switches
			Switches="--no-gui --eval"
			;;
    esac
fi

## Make sure MATLAB/Octave know the intent
# note: the export is required
export CONTINUOUS_INTEGRATION=true
export CI=true

## Actually run the test suite
cd unitTest
TESTDIR=`pwd`
# also CD in MATLAB/Octave to make sure that startup files
# cannot play any role in setting the path
${Runner} ${Switches} "cd('${TESTDIR}'); runMatRadTests"
exitIfError $?
cd ..

## Post-processing

# convert MD report into HTML using pandoc if available
MDFILE="unitTest/results.test.md"
if [ ! -z `which pandoc` ]; then
	if [ -f $MDFILE ]; then
	    HTMLFILE=${MDFILE/md/html}
	    # replace the emoji while we're at it
	    pandoc -f markdown -t html $MDFILE -o $HTMLFILE
	    sed -i -- 's/:heavy_exclamation_mark:/❗️/g' $HTMLFILE
	    sed -i -- 's/:white_check_mark:/✅/g' $HTMLFILE
	    sed -i -- 's/:grey_question:/❔/g' $HTMLFILE
	fi
fi

