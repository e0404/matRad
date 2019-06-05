echo off
REM Make sure some failures are detected by the CI runners
REM function exitIfError {
	REM # pass "$?" as argument: i.e. the exit status of the last call
	REM if [ "$1" -ne 0 ]; then
		REM exit $1;
	REM fi
REM }

REM Handle Runner variables
set Runner=%1

REM if [ -z "$Switches" ] ; then
	REM case "$Runner" in
		REM *matlab* )
			REM Switches="-nodesktop -r"
			REM ;;

		REM *octave* )
			REM Switches="--no-gui --eval"
			REM ;;

		REM * )
			REM # Fall back to Octave switches
			REM Switches="--no-gui --eval"
			REM ;;
    REM esac
REM fi

REM Make sure MATLAB/Octave know the intent
REM note: the export is required
set CONTINUOUS_INTEGRATION=true
set CI=true

REM Actually run the test suite
cd unitTest
set TESTDIR="%cd%"
REM also CD in MATLAB/Octave to make sure that startup files
REM cannot play any role in setting the path
%Runner% --eval "cd('%TESTDIR%'); matRad_runTests" > ../runtests.log
REM exitIfError $?
cd ..

