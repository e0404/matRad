# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Userfolders can now also be set via environment variable `MATRAD_USERDATA`
- Documentation: Documented the userfolder feature and its usage as well as other datastructures more clearly

### Fixed
- possible negative doses in finesampling engine due to extrapolation in kernel interpolation
- correct parsing of all optional arguments of the `traceCube` function for `matRad_RayTracer`

## [3.2.2]

### Fixed
- Fixed `matRad_version` in case of tagged releases.
- Fixed documentation: small documentation correction.

## [3.2.1]

This patch fixes a multitude of issues with the new MATLAB Desktop from R2025 and reported minor issues in import/export, helper functions, and other utility functions. Apart from that, it introduces new flexibilities under the hood, such as single precision and raytracer vectorization, that do not intentionally change outward-facing behavior and can be tested before broader behavioral changes in a future major release. The license was also changed to BSD 3-Clause.

### Added
- Dose engines can now optionally run calculations in single precision while the default remains double. A new `precision` configuration property controls this.
- GPU acceleration is available as an opt-in property for optimization. Helper functions for translating matRad data structures to and from GPU arrays were added.
- The Siddon raytracer is now implemented as a class with vectorized ray processing and optional single-precision forcing.
- Pencil beam engines now expose a `traceOnDoseGrid` switch, defaulting to `false`, to optionally retain radiological depth cubes on the dose grid.
- Improvements to the phantom builder now also allow mm coordinates for phantom definition.
- The FRED interface was updated with new test data, improved version compatibility, and the ability to force `ijFormatVersion`.
- DICOM import now imports passively scattered proton beams, including gantry and couch angles.
- Documentation: added the full Sphinx and Read the Docs documentation build pipeline, including `readthedocs.yml`, a GitHub Actions workflow, the `docs/` folder, and the reStructuredText documentation structure.
- CI: added a standalone build step to GitHub Actions workflows with a matrix build for Windows, Linux, macOS Intel, and macOS Apple Silicon.
- CI: added preliminary pre-commit hook configuration with `miss_hit` and `codespell`, not yet enforced.
- CI: added a GitHub Actions workflow for documentation building triggered by changes to `docs/`.

### Changed
- Streamlined sequencing and 3D conformal calculations.
- Optimizer instantiation was reworked to allow more configuration options via `propOpt`.
- The `finalizeDose` call in dose engines was moved to `calcDoseForward` and `calcDoseInfluence`.
- The matRad license changed to BSD 3-Clause.
- Documentation: updated docstrings across many files to be Sphinx Napoleon compatible.
- Documentation: updated copyright notices to 2026.
- CI: made the coverage report workflow more tolerant to errors.
- CI: updated the MOcov submodule to include an md5 fix.

### Fixed
- Variance calculation from MC statistics can now be computed correctly.
- `matRad_plotSlice` input parsing was improved, including a fix for empty figure opening due to colormap array requests.
- TOPAS now correctly supports multiple alpha/beta values.
- Fixed range-shifter handling issues in MC dose calculation interfaces.
- Fixed a typo in the RBE model fallback load path.
- Fixed a typo in `addMUdataFromMachine`.
- Corrected the DICOM attribute for `SliceLocation`.
- Fixed a slight dimension interpretation issue in `cubeIndex2worldCoords`.
- Fixed scenario listing and the robustness field when serializing objectives to structs or displaying them in the CST.
- Multiple GUI fixes were added for MATLAB 2025 compatibility.
- GUI fixes include a missing plot handle, empty figure handles returned when the GUI is globally disabled, a `plotSlice` colormap issue, and scrolling in the viewing widget under Octave when `CurrentPoint` is empty.
- `numOfbeams` is no longer required because it can be inferred.

## [3.2.0]

### Added
- Added the FRED MC interface, if installed.
- Added VHEE planning with a generic unfocused beam and a focused beam. The generic beam can also be forwarded to TOPAS.
- Added a new `matRad_plotSlice` function with keyword/value syntax for more intuitive slice plotting.
- Added new examples for usage of FRED and VHEE and a workflow example for comparing dose calculation on synthetic CT to planning CT.
- CI: added a new `.gitlab-ci.yml` file to support GitLab CI/CD, including test and package stages, artifact handling, and configuration for MATLAB container images and licensing.
- CI: added more comprehensive dose calculation tests.
- Project: added a `.gitattributes` file to standardize line endings, treat certain file types as binary, and ensure `.m` files are not marked as executable.
- Project: added new contributors.

### Changed
- Updated examples to use `matRad_plotSlice`.
- The analytical functions from the Bortfeld Bragg Peak Model are now public and can be used to compute standard approximations such as range-energy relationships.
- CI: added `Global_Optimization_Toolbox` to the MATLAB products list in `.github/actions/test-matlab/action.yml`.
- CI: made the coverage PR comment step in `.github/workflows/coverage-report.yml` tolerant to errors to avoid workflow failures.

### Fixed
- The DICOM import widget now allows selection of multiple RTDose files.
- The DICOM import widget and importer now handle selected patients more consistently and robustly.
- The DICOM exporter writes quantities beyond dose, and the importer now tries to import them correctly.
- The DICOM exporter now always writes `ReferencedRTPlanSequence`, and the importer can now survive without it.
- The DVH widget no longer throws a warning during updates and handles scenarios and missing `xlabel` `axesHandle` parameters more robustly.
- GUI fixes were made for setting gantry angles and other parameters in the PlanningWidget.
- `EXTERNAL` contours are now recognized correctly.
- Improved performance when obtaining the Jacobian structure in optimization.
- Available classes such as dose engines are now cached for faster loading.
- GUI fixes were added for use in MATLAB Online.

## [3.1.0]

### Added
- Introduced a major file structure overhaul into organized subfolders such as `matRad`, `thirdParty`, and `examples` to improve clarity and maintainability.
- Introduced the `userdata` folder to maintain custom data.
- Introduced comprehensive scenario management, including support for 4D phase scenarios and automated scenario model instance tests.
- Added multiple robust optimization methods, including COWC, OWC, VWWC, and expected value.
- Transitioned dose calculation to object-oriented dose engines, improving structure and maintainability.
- Added a customizable TOPAS interface for ions and an experimental interface for photons.
- Reworked the workflow of existing Monte Carlo interfaces into the new engine format.
- Added a new coordinate system handling approach with separate world and cube systems and dedicated transformation functions.
- Recalculated the `proton_Generic` machine with stored phase space parameterization to facilitate consistent MC and PB calculations.
- Added a generic helium dataset including LET.
- Added a LET-based helium model.
- Added multiple variable RBE models for protons and helium with an object-oriented data model architecture.
- Added BED optimization.
- Replaced MATLAB GUIDE with a modern widget-based GUI.
- Added light and dark mode.
- Added a DICOM exporter for CTs, RTStruct, safeguarded RTPlan for photons, and RTDose.
- Refactored the DICOM importer for better use from scripts.
- Introduced MOxUnit and MOcov for automated unit tests, including example tests, as submodules.

### Changed
- Large parts of the GUI are now Octave compatible.
- Default configuration options are now stored in `MatRad_Config` under the `defaults` struct, with a compatibility layer for older access patterns.
- `matRad_calcCubes` was changed to accept a variety of Monte Carlo-related fields without changing typical current usage.
- Biological models are now defined in a fundamentally different way, and downward compatibility is not guaranteed.
- The object-oriented scenario models and biological models may cause issues in older scripts if matRad cannot infer the models.
- The old dose calculation functions remain in a compatibility and deprecation layer, but some configuration options may not work as intended.
- Unit tests are now run as GitHub Actions on MATLAB R2022b, the latest MATLAB release, and Octave 6.
- Compatibility is tested for Octave 6 through 9.
- Improved performance in interpolation.
- Improved performance in optimization.
- Cleaned up code for more consistent use of `MatRad_Config` error, warning, and logging mechanisms.
- Octave compatibility is tested for Octave 6 through 9, although IPOPT still needs to be compiled individually.
- The isocenter is now always given in world coordinates corresponding to the CT plane coordinates, whereas it previously used its own cube coordinate system based on voxel index times resolution.
- Other coordinate system bug fixes may induce changes when an existing script is rerun.
- Starting from major version `3.*`, matRad follows semantic versioning in `major.minor.patch` form.
- Major releases include major new features, such as a new modality, and may break top-level API compatibility.
- Minor releases include smaller new features, such as a new optimizer, objective, biomodel, or dose calculation algorithm, while preserving compatibility within the major release.
- Patch releases fix bugs and do not generally introduce new features, except for minimal configuration options that mitigate bugs in special cases.

### Fixed
- Resolved issues with `ompMC` mex file compilation and Octave compatibility warnings.
- Corrected path issues and file handling, especially for temporary directories and submodules.
- Fixed bugs in optimization objectives and constraints for special input cases.
- Fixed issues in DICOM import that expected non-standard tags.

### Deprecated
- The previous procedural dose calculation workflow has been superseded by object-oriented dose engines, although compatibility layers remain in place.

### Removed
- MATLAB GUIDE-based GUI usage as the primary GUI architecture.

## [2.10.1]

Release with small updates, cleanups, and bug fixes.

### Added
- Added blue/white/red difference maps to the available colormaps.
- Added the option `pln.propDoseCalc.useGivenEqDensityCube`, defaulting to `false`, to directly use the literal values from `ct.cube` and omit HU to WEQ conversion from `ct.cubeHU`.
- Added the option `pln.propDoseCalc.ignoreOutsideDensities`, defaulting to `true`, to disable or enable inclusion of WEPL outside the patient contour during ray tracing.

### Changed
- Removed the hardcoded penumbra width in photon dose calculation so it can now be stored in the machine file as `machine.data.penumbraFWHMatIso`.
- Updated `ompMC` to use a virtual Gaussian source with measured penumbra values, including precompiled mex files.
- Updated Travis CI testing, speeding it up by using precompiled mex interfaces and including testing with MATLAB on Ubuntu, with Azure DevOps as fallback.
- Added GitHub automation such as Stalebot and issue and PR templates.
- Optimized Jacobian evaluation with a reported 10x to 100x speedup.
- Documentation: updated code documentation.

### Fixed
- Fixed a 3D view issue caused by inconsistent angles in `pln` and `stf`.
- Fixed incorrect DICOM UIDs and writing order in DICOM export.
- Fixed a colormap issue in plotting.
- Added new handling of environment checking with `matRad_cfg`, while retaining the old function.
- Removed unnecessary `global` statements before `matRad_cfg`.

## [2.10.0]

Second release of matRad. Despite major incompatibilities with "Alan", the project kept the major version number `2` to preserve a more consistent versioning scheme going forward. The team thanks all new contributing authors listed in `AUTHORS.txt`.

### Added
- Added integration tests using Travis CI with Octave, excluding GUI functionality.
- Added the `matRad_rc` script to configure matRad paths.
- Added version printing through `matRad_version`, and the version is now shown in the GUI and when using `matRad_rc`.
- Added seven new MATLAB example scripts demonstrating matRad usage.
- Added basic interfaces to the open-source photon and proton MC engines `ompMC` and `MCsquare`.
- Overhauled the optimization interface using object-oriented programming and integrated the `fmincon` optimizer from the MATLAB Optimization Toolbox.
- Changed the `cst` variable and added a conversion script for old to new `cst` structures in the `tools` folder.
- Separated CT, optimization, and dose calculation grids to allow different resolutions.
- Added developer mode and simplified educational mode for the graphical user interface.
- Organized base data and default phantoms into subfolders.
- Added DICOM export. The `dicomImport` folder was renamed to `dicom`.
- Added a DICOM import GUI.
- Added binary import and export functionality in both script and GUI form.
- Overhauled the standalone project file.
- Added a standalone toolbox file for matRad.
- Dose calculation now uses generalized initialization scripts.
- Added `matRad_compareDose` to compare dose distributions with difference and gamma analysis.
- Added more tools for visualization and data analysis in the `tools` folder.
- Added support for defining a range shifter.
- Added quality indicator and DVH display wrapper scripts.
- Added a wrapper to allow 3D conformal planning using `dij` collapsing.
- Added new colormap handling to allow integration of custom colormaps.
- Modularized slice display through dedicated plotting functions, including generation of 3D views.
- Added the global configuration object `matRad_cfg` backed by `MatRad_Config` to store default values and provide a logging interface.

### Fixed
- Many bug fixes and many new bugs.

## [2.1.0]

First official release of matRad.

### Added
- Added the IPOPT optimizer for constrained optimization.
- Added validated ray tracing.
- Added validated pencil beam particle dose calculation.
- Added validated singular value decomposed pencil beam photon dose calculation.
- Added DICOM import including dose and particle pencil beam scanning plan objects.
- Added a standalone matRad version.
- Improved the GUI workflow.

### Fixed
- Many bug fixes and many new bugs.
