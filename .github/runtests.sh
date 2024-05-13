#!/usr/bin/env bash


## Make sure some failures are detected by the CI runners
function exitIfError {
	# pass "$?" as argument: i.e. the exit status of the last call
	# currently octave 6 can finish with a segfault when the program is closed due to some bug, for now we try to ignore it
	if [ "$1" -ne 0 ] && [ "$1" -ne 139 ]; then
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
			Switches="-batch"
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
# also CD in MATLAB/Octave to make sure that startup files
# cannot play any role in setting the path
${Runner} ${Switches} "back=cd('MOxUnit/MOxUnit'); moxunit_set_path; cd(back); cd('matRad'); matRad_rc; moxunit_runtests('test','-recursive','-junit_xml_file','testresults.xml'); exit(double(~ans))" > ../runtests.log #2> ../runtests.err put stdout to log, but let it print error messages
exitIfError $?

