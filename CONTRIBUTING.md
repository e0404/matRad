You want to contribute to matRad? You want to submit a bug report, improve the documentation, have feature requests or write code that can be incorporated into matRad? Great! matRad is a constantly evolving software that is waiting for people like you who want to improve it. 

Please read and follow the subsequent guidelines. This will help us to address your issue, assess your changes and finalize your pull request as efficiently as possible. 

## Ground rules

### Before reporting a bug or asking questions 

If you encounter problems with matRad, please consider the following guidelines before submitting issues on our github page. 

- Check you are using the newest version of matRad.
- Check the description of how to set up matRad and its technical documentation in the [wiki](https://github.com/e0404/matRad/wiki). 
- Go through the relevant workflow examples and see if they answer your question. The workflow examples can be found in the subfolder [_examples_](https://github.com/e0404/matRad/tree/master/examples). 
- Check open and closed issues for your question.

It is very helpful for us, if you search for the error by yourself a bit before asking. You might be able to fix the error by yourself. Please keep in mind that we are no support team but researcher. 

### How to report a bug and ask questions

Still having problems? Use the GitHub functionalities to file a new issue, or directly fork the matRad repository and send pull requests with your own custom developments. Please take into account the given rules given below and be patient! 

- Provide the following information:
	- OS
	- Software environment and version (executing `[env, versionString] = matRad_getEnvironment` tells you the software environment matRad is running on) 
	- a **minimum example** of your attempted workflow / what causes the problems
- Post only questions related to matRad.
- Don’t write an email if you want to report a bug. Write only an email if you want to cooperate or plan a bigger project with matRad. 

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
	
### matRad programming style

If you add or change code, follow the matRad programming style:

**White space**

Always use white spaces before and after the "=" sign upon assigning new values to variables, e.g.

`thisIsA = validAssignment;`

**Variable names**

matRad variable names start with a lower case letter. Concatenated variable names indicate a new word by an upper case letter e.g.

`newVariableName = 42;`

**Absolute paths**

At some points there are issues with relative paths to matRad files which are annoying (e.g. within the unit testing). The problems occur when you run certain matRad functionality and are currently not in the matRad root directory. To avoid this, we should make sure, that

- all references to files are absolute paths, which are created at runtime through locating the file that is executing the code (e.g. via `mfilename('fullpath')`). Alternatively, `which` could be used.
- if something needs to be executed within a directory, also use the absolute path to `cd` into the directory. But before, use `pwd` to get the users directory and after the execution `cd` back into the directory.
	
## Cooperating or planning a bigger project with matRad

If you want to cooperate or plan a bigger project with matRad write us an [email](mailto:contact@matRad.org).
