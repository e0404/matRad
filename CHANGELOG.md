# Changelog

## Development Changes

### Bug Fixes
- DICOM Import widget allow selection of multiple RTDose files.
- DICOM Import Widget and importer handle selected patient more consistently and robustly. 
- DICOM Exporter writes quantities beyond dose, importer tries to import them correctly.
- DICOM Exporter now always writes ReferencedRTPlanSequence. Importer can now survive without it.
- DVH widget does not throw a warning in updates, handle scenarios correctly / more robustly and missing xlabel axesHandle parameter.


## Version 3.1.0 - "Cleve"

### Major Changes and New Features

#### File Structure Overhaul 
- Major restructuring of files into organized subfolders, such as matRad (core code), thirdParty, examples, etc., to improve clarity and maintainability.
- Introduction of userdata folder to maintain custom data

#### Scenario Management and robust / 4D optimization
- Introduced comprehensive scenario management (scenario models), including support for 4D phase scenarios and automated scenario model instance tests.
- Multiple robust optimization methods (COWC, OWC, VWWC, expected value) 

#### Object-Oriented DoseEngines & new Monte Carlo interfaces 
- Transitioned from procedural dose calculation to an object-oriented approach, significantly improving the structure and maintainability of the dose engines.
- Added customizable TOPAS interface for ions (and experimental for photons)
- Workflow of the existing Monte Carlo interfaces has been completely overhauled in the new engine format
- New handling of coordinate system (separation into world / cube systems with dedicated transformation functions) to ease readability
- Recalculated proton_Generic machine with stored phase space parameterization to facilitate consistent MC / PB calculations.

#### Helium planning
- matRad now contains a Generic helium dataset including LET
- LET-based Helium model

#### Extended biological modeling
- Multiple variable RBE models for protons and helium added with an Object-oriented Datamodel Architecture
- BED optimization

#### Widget-Based GUI
- Replaced Matlab's GUIDE-based approach with a modern widget-based GUI.
- Large parts of the GUI are now Octave Compatible
- Light & Dark Mode

#### DICOM Exporter and Importer
- added a DICOM exporter for CTs, RTStruct, RTPlan (photons, safeguarded) and RTDose
- Refactored the DICOM importer for better use from scripts

#### Possibly (and probably) Breaking Changes to matRad core workflow and functions
While we try to keep downwards compatibility (and will provide fixes if breaking changes are detected), here are some potential/probable dealbreakers
- The coordinate handling of the isocenter changed. The isocenter is now always given in "world" coordinates (i.e., corresponding to the ct plane coordintes). Before, the isocenter resided in its own "cube" coordinate system (voxel index * resolution)
- Some other coordinate system bug-fixes might induce changes when an existing script is rerun
- Default configuration options now stored in MatRad_Config under "defaults" struct. There is a compatibility layer, but this might break under user changes
- Changed matRad_calcCubes to accept a variety of different fields for Monte Carlo, without changing the current usage
- Biological Models are defined in a completely different way now and downwards compatibility is not guaranteed.
- The object oriented scenario models and biological models could procude issues in old scripts if matRad can not infer the models
- While the old dose calculation functions have been kept in a compatibility / deprecation layer, some configuration options might not work as intended

### Other Enhancements, Documentation, and Testing

#### New unit-testing framework
- Introduced Usage of MOxUnit and MOcov for automated unit tests (and the example tests). They are included as submodules.
- Unit tests now runnig as GitHub Actions on Matlab R2022b, the latest release, and Octave 6

#### Improved Octave Compatibility:
- Compatibility tested for Octave 6 to 9. 
- Octave compatibility not always optimal, and IPOPT needs to be compiled individually.

#### Performance Improvements & Code Cleanup
- Performance improvements and updates on interpolation
- Performance improvements in optimization
- Code cleanup for consistent use of MatRad_Config's error / warning / logging mechanism.

### Bug Fixes
- *Compilation Fixes:* Resolved issues with ompMC mex file compilation and Octave compatibility warnings.
- Corrected path issues and file handling, especially for temporary directories and submodules.
- Fixed some bugs in optimization objectives & constraints for special input cases
- Fixed issues in DICOM import expecting non-standard tags

### Semantic Versioning
- Starting from major version 3.*, matRad will follow rigorous semantic versioning in a major.minor.patch style
  - **Major:** Major releases include major new features (e.g. a new modality) and/or do not guarantee downwards compatibility of the top-level API, which we consider calls to all functions on the top-level of the "matRad" folder. Since a lot of configuration of these functions can be done per `pln.prop*` and other propertie sin `pln`, this might happen more quickly than one might think.
  - **Minor:** Minor releases incldue minor new features (e.g. a new optimizer, objectives, biomodel or dose calculation algorithm). Downwards compatibility (within the major release) is preserved.
  - **Patch:** Patch versions only fix bugs and do not introduce new features. Exceptions could be the exposure of new, minimal configuration options to mitigate a bug occuring in special circumstances.

## Version 2.10.1 - Patch release for "Blaise" 
Release with small updates, clean-ups and bugfixes    
- Bugfix in 3D view due to inconsistent angles in pln & stf
- Bugfix for using incorrect dicom UID's and wrong writing order in the dicom export
- Bugfix for weird colormap issue in plotting
- New handling of environment checking with matRad_cfg (old function is still working)
- Code documentation update
- Remove hardcoded penumbra width in photon dose calculation -> can now be stored in machine file (machine.data.penumbraFWHMatIso)
- Update to ompMC to use virtual Gaussian source (uses measured penumbra value) incld precompiled mex files
- Remove useless global statements before matRad_cfg
- Add blue/white/red difference map to colormaps (in the correct way)
- Updated TravisCI testing (Sped up by using pre-compiled mex interfaces and including testing with Matlab (on Ubuntu), Azure DevOps as fallback
- Github gimmicks added: Stalebot, Issue & PR Templates
- Code optimization for jacobian evaluation (x10-100 speedup)
- New option pln.propDoseCalc.useGivenEqDensityCube (default false) to directly use the literal values from ct.cube and omit HU to WEQ conversion from ct.cubeHU
- New option pln.propDoseCalc.ignoreOutsideDensities (default true) to disable/enable inclusion of WEPL outside the patient contour in ray-tracing

## Version 2.10.0 - "Blaise"
Second Release of matRad. Note that despite major incompatibilities with "Alan" we still chose major version number "2" to have a consistent versioning in the future.
We want to thank all new contributing authors (see AUTHORS.txt)
The new release contains:
- Integration tests using TravisCI (with Octave, so no GUI functionalities are tested)
- matRad_rc script to configure matRad paths
- matRad version can now be printed with matRad_version, version correctly shown in GUI and when using matRad_rc
- Seven new Matlab example scripts to demonstrate use of matRad 
- Added basic interfaces to the open-source photon/proton MC engines ompMC/MCsquare
- Overhaul of the optimization interface using OOP and integration of the fmincon optimizer from Mathworks' MATLAB Optimization toolbox
- Changes to the cst variable (new script to convert old to new cst structures in tools folder)
- Separation of ct/optimization/dose calculation grids to allow different resolutions
- The graphical user interface can now be started in developer mode (more error information) and a simplified educational mode (some functionalities disabled, less buttons)
- Base data and default phantoms now organized in subfolders
- DICOM export added (only script, dicomImport folder renamed to dicom)
- DICOM import GUI
- Binary import and export functionalities (script & GUI)
- Overhauled the standalone project file
- Standalone toolbox file for matRad
- Dose calculation now uses generalized initialization scripts
- matRad_compareDose tool to compare two dose distributions with difference and gamma analysis
- More tools for visualization and data analysis in the tools folder
- Possibility to define range shifter
- Quality indicator & DVH display wrapper scripts
- Wrapper to allow 3D conformal planning using dij collapsing
- New colormap handling to allow integration of custom colormaps
- Modularization of slice display by dedicated functions in plotting folder including generation of 3D views
- New global configuration object (matRad_cfg <- MatRad_Config.m) to store default values and with logging interface
- Many bug fixes and many new bugs..	

## Version 2.1 "Alan"
First official release of matRad including
- New optimizer IPOPT for constrained optimization
- Validated ray tracing
- Validated pencil beam particle dose calculation
- Validated singular value decomposed pencil beam photon dose calculation
- DICOM import including dose and particle pencil beam scanning plan objects
- matRad standalone version
- Improved GUI workflow
- Many bug fixes and many new bugs...
