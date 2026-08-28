You want to contribute to matRad? You want to submit a bug report, improve the documentation, have feature requests or write code that can be incorporated into matRad? Great! matRad is a constantly evolving software that is waiting for people like you who want to improve it.

Please read and follow the subsequent guidelines. This will help us to address your issue, assess your changes and finalize your pull request as efficiently as possible.

## Ground rules

### Before reporting a bug or asking questions

If you encounter problems with matRad, please consider the following guidelines before submitting issues on our github page.

- Check you are using the newest version of matRad.
- Check the description of how to set up matRad and its technical documentation in the [wiki](https://github.com/e0404/matRad/wiki).
- Go through the relevant workflow examples and see if they answer your question. The workflow examples can be found in the subfolder [_examples_](https://github.com/e0404/matRad/tree/master/examples).
- Check open and closed issues for your question.

It is very helpful for us, if you search for the error by yourself a bit before asking. You might be able to fix the error by yourself. Please keep in mind that we are no support team but researchers that can only point you in the right direction. The biggest joy when using research software lies in figuring out details yourself, which also gives you the biggest learning effect and even might enable you to contribute something to matRad which can the be used by other researchers in the future.

### How to report a bug and ask questions

Still having problems? Use the GitHub functionalities to file a new issue, or directly fork the matRad repository and send pull requests with your own custom developments or bug fixes. Take into account the given rules given below and be patient!

- Provide the following information:
	- OS
	- Software environment and version (executing `[env, versionString] = matRad_getEnvironment` tells you the software environment matRad is running on). Also, when running `matRad_rc` (since version 2.10.0), you get full printed information about version (and branch/commit if you are working with git).
	- a **minimum example** of your attempted workflow / what causes the problems. Did you cchange any matRad files?
- Post only questions related to matRad.
- Please don’t write an email if you want to report a bug. An issue on github is more helpful to others and will probably also speed up the process of fixing the bug.

### Adding functionalities or changing code

If you want to add functionalities or change code:

1. Create your own fork of the code.
2. Do the changes in your fork.
3. In case you like the changes and think they should be available for everyone:
	- Be sure you followed the matRad programming style.
	- Search GitHub for an open or closed Pull Request that relates to your submission in order to avoid duplicate effort.
	- Commit your changes using a descriptive commit message.
	- Push your branch to GitHub.
	- In GitHub, send a pull request to the [dev branch](https://github.com/e0404/matRad/tree/dev).

A more detailed description of the repository layout, architecture and development workflow (also aimed at AI coding agents) can be found in [AGENTS.md](AGENTS.md).

### Contribution checklist

- **Tests:** matRad uses [MOxUnit](https://github.com/MOxUnit/MOxUnit); run `matRad_runTests` (see [test/README.md](test/README.md)). Please add tests for new functionality and a regression test for every bug fix, preferably before writing the fix. Tests need to pass on the stable MATLAB release (currently R2022b), the latest MATLAB release and Octave (currently version 8), which are all checked in CI.
- **Octave compatibility:** keep code working in Octave where possible. Compatibility wrappers live in `matRad/util/octaveCompat/`.
- **Object-oriented design:** new algorithms should be implemented as classes in the existing hierarchies (dose engines, stf generators, sequencers, optimizers, scenario and biological models, ...) and configured through `pln.prop*`; the procedural top-level `matRad_*` functions only wire these objects together.
- **Code style:** style, naming and complexity rules are defined in `.miss_hit` and checked with [MISS_HIT](https://florianschanda.github.io/miss_hit/) via [pre-commit](https://pre-commit.com/) (`pip install pre-commit miss_hit`, `pre-commit install`). The hooks check complete files, so touching legacy code can trigger extensive reformatting. You may skip the hooks (`git commit --no-verify`) to keep your pull request focused, but the code you add or change must satisfy the rules.
- **Changelog:** note every user-relevant change in `CHANGELOG.md` under *Unreleased*, following [Keep a Changelog](https://keepachangelog.com/) conventions (the file is compatible with [kacl](https://github.com/mschmieder/python-kacl)).
- **Documentation:** update the Sphinx documentation in `docs/` (build with `pip install -r docs/requirements.txt` and `sphinx-build docs _doc`) when behavior, options or data structures change. The API reference is generated from the docstrings.
- **Authorship:** if you contribute code for the first time, add your name to `AUTHORS.txt` and to the author list in `CITATION.cff`.

### matRad programming style


If you add or change code, follow the matRad programming style:

**White space**

Always use white spaces before and after the "=" sign upon assigning new values to variables, e.g.

`thisIsA = validAssignment;`

**Variable names**

matRad variable names start with a lower case letter. Concatenated variable names indicate a new word by an upper case letter e.g.

`newVariableName = 42;`

**Console output and default parameters with MatRad_Config**

Since matRad 2.10.0, it uses a configuration class `MatRad_Config` to enable setting default parameters and to control the console output. The configuration class uses the "Singleton"-programming paradigm, i.e. only one instance of it exists at any time so that the default parameters and the logging mechanism stay consistent while running matRad.
If you write information to the console, use the respective functions within the class. The configuration object can be instantiated at any point in a function by calling ```matRad_cfg = MatRad_Config.instance();```

**Absolute paths**

At some points there are issues with relative paths to matRad files which are annoying (e.g., within the unit testing). The problems occur when you run certain matRad functionality and are currently not in the matRad root directory. To avoid this, we should make sure, that

- all references to files are absolute paths, which are created at runtime through locating the file that is executing the code (e.g. via `mfilename('fullpath')`). Alternatively, `which` could be used.
- if something needs to be executed within a directory, also use the absolute path to `cd` into the directory. But before, use `pwd` to get the users directory and after the execution `cd` back into the directory.


## You want to cooperate with us on some research project?

If you plan a research project that uses matRad and you would like to cooperate do not hesitate to write us an [email](mailto:contact@matRad.org).
