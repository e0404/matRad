## matRad Tests
matRad tests run with [MOxUnit](https://github.com/MOxUnit/MOxUnit). Coverage is generated with [MOcov](https://github.com/MOcov/MOcov).

matRad runs through (most of the) examples (each one counting as one test) on top of defined unitTests using MOxUnits Syntax.
Multiple Matlab versions as well as Octave are tested, all running on Ubuntu.

### Testing Locally

1. Clone MOxUnit into a *separate* folder (not within or in subfolders of the matRad root dir).
2. Navigate into the MOxUnit code folder (within the MOxUnit) and call `moxunit_set_path`.
3. Navigate back into the matRad root folder. Make sure your matRad folder is on the path by running `matRad_rc`.
4. Run tests on matRads test folder by calling `moxunit_runtests('test','-recursive');`. Test output will be written to the command line. 
    - If you want to create a html-formatted output, you can also create a logfile by calling `moxunit_runtests('test','-recursive','-logfile', 'tests.log')`
    - You can limit your tests to a subfolder only (e.g., if you wrote new unit tests and want to do isolated bugfixing) by running on the subfolder directly, for example: `moxunit_runtests('test/tools','-recursive');`