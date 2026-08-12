# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added elastic image registration for 4D CT scenarios: `matRad_DemonsImageRegistration` (based on `imregdemons`, requires the Image Processing Toolbox) computes pull or push deformation vector fields between a reference CT scenario and the remaining scenarios, stored in `ct.dvf` / `ct.dvfMetadata` in mm by default and compatible with `matRad_doseAcc`, and can propagate contours to the other scenarios using push fields. New registration algorithms can be added by subclassing `matRad_ImageRegistrationBase`. Also adds `matRad_checkEnvImageProcessingRequirements` to check availability of the Image Processing Toolbox (MATLAB) / image package (Octave)
- Added helper to create ring VOIs from margins around existing structures: `matRad_createRing`, optionally grown around the union of the base VOI over all CT scenarios (e.g. a ring around the full motion of a 4D target)
- `matRad_PhantomVOISphere` accepts an optional `innerRadius` to create spherical shells, exposed via the new `matRad_PhantomBuilder.addSphericalShellOAR`
- Support for VMAT (volumetric modulated arc therapy) planning of photon plans, enabled with `pln.propOpt.runVMAT`. The `PhotonVMAT` stf generator places a hierarchy of fluence-optimization angles, DAO control points and dose-calculation angles along one or more arcs (`pln.propStf.maxFMOGantryAngleSpacing` / `maxDAOGantryAngleSpacing` / `maxGantryAngleSpacing`, `minAperturesPerFMOBeam`, `arcIndex`); fluence optimization runs on the FMO angle subset; the siochi sequencer decomposes each fluence map and spreads the resulting apertures over that beam's DAO control points; and direct aperture optimization then optimizes aperture weights, leaf positions and gantry times subject to the machine's gantry rotation speed, leaf speed and monitor unit rate limits (`machine.constraints`), optionally with interpolated leaf positions between control points (`pln.propSeq.continuousAperture`) and Jacobi preconditioning (`pln.propSeq.preconditioner`). `matRad_calcDeliveryMetrics` reports the delivery of the resulting plan and `matRad_refineApertureArc` recalculates the dose on a finer gantry angle grid. Demonstrated in the new `matRad_example22_photonsVMAT` and described in the new *Sequencing, Direct Aperture Optimization & VMAT* section of the technical guide
- Fluence objectives, a family of optimization functions acting on the bixel weight vector directly instead of on a dose related quantity, so that they take no part in the backprojection. `matRad_OptimizationFunction` is the new common ancestor of the existing `matRad_DoseOptimizationFunction` (unchanged in name and API) and of the new `matRad_FluenceOptimizationFunction`. `FluenceObjectives.matRad_FluenceVariance` penalizes the fluence variance within each beam, which limits the modulation a plan may use and is what makes an arc sequenceable into few apertures; it needs only the bixel-to-beam mapping from the dij, so it works for photons and particles alike. Fluence objectives are given either in `pln.propOpt.fluenceObjectives` or alongside the dose objectives in `cst{i,6}`; direct aperture optimization warns and ignores them, since smoothing a sum of deliverable apertures has no well defined meaning
- VMC++ photon Monte Carlo is available as a regular dose engine (`matRad_PhotonVmcEngine`, short name `VMC`) instead of a standalone script, so it is selected and configured like every other engine through `pln.propDoseCalc`
- `matRad_visApertureInfo` gained a view-first API (`matRad_visApertureInfo(apertureInfo, view, ...)`) and additional views: `perBeam` (tabbed, one aperture grid per beam), `trajectory` (leaf trajectories over the arc, with leaf pairs labelled by index or drawn at their physical coordinates), `interactive` and `animate` (arc playback in delivery time, with the leaves and gantry moving between control points)

### Fixed

- Progress display in the console no longer gets mangled when other log output is written in between updates: `matRad_progress` accepts an optional `linereset` argument, used by the pencil-beam dose engines to restart the progress display for each beam
- Removed a duplicated progress printout per beam during stf generation
- The debug visualization of external beam stf generation (`visMode > 0`) crashed on undefined variables since the class-based generator refactor and now works again; the LPS subplot also gets its own axis labels instead of relabeling the beam's eye view
- Direct aperture optimization: the constraint Jacobian was wrong for several entries and the VMAT leaf position bookkeeping was inconsistent between the aperture vector and the aperture info struct; both are now covered by finite-difference consistency tests (`test_daoGradientConsistency`)
- `matRad_recalcApertureInfo` populated `nextDAOBeamIx` from `lastDAOBeamIx`, so every recalculated beam believed its next DAO anchor was its previous one and the interpolation between control points was weighted against the wrong pair of apertures
- VMAT arc pruning and refinement produced invalid geometry for arcs traversed in decreasing gantry angle direction, and the refined apertures were written to the wrong shape fields
- `matRad_calcMaxLeafSpeed` computed the right-leaf speed from the left-leaf positions, so right-leaf motion never contributed to the reported maximum leaf speed
- VMC++ dose calculation did not run on Windows and used inconsistent coordinate conventions
- The MCsquare Windows binary is called with an explicit path prefix, so it is found regardless of the current directory
- `MatRad_Config.loadobj` merged saved structs using the wrong loop variable
- The Sphinx documentation now builds without warnings: the parameter tables of several docstrings were malformed reStructuredText and rendered as run-together prose, five functions had a comment banner in place of their summary line, and three `.rst` files had a broken cross reference, a short title underline and a misindented list
- Matlab's two-argument `round(x,n)` was used in VMAT steering information generation, the sampling report and the DVH band plot, which errors on Octave because its `round` takes the value alone. These now call the new `matRad_roundCompat`, and a test guards against reintroducing the two-argument form
- Classes living in `@` folders were missing from the API documentation and the MATLAB module index, since `automodule` does not descend into them: `matRad_OptimizationProblem` and its DAO and VMAT subclasses, `matRad_DoseEngineBase` and `matRad_ParticleFREDEngine` are now documented explicitly. The `FluenceObjectives` package was likewise missing and is now included

### Changed

- Consistency handling of gantry and couch angles in external beam stf generators was reworked: the property setters no longer silently pad or trim the respective other angle vector. Inconsistent numbers of gantry and couch angles now throw an error, both when set via `pln.propStf` and (on use) when set directly on the generator. As a documented convenience, a scalar couch angle is valid for any number of gantry angles and is applied to all beams
- Refactored leaf/spot sequencing into an object-oriented class hierarchy (`matRad_SequencerBase` with the photon MLC sequencers `matRad_PhotonLeafSequencerSiochi`/`matRad_PhotonLeafSequencerXia`/`matRad_PhotonLeafSequencerEngel` and the particle spot sequencer `matRad_ParticleScanningSequencerSpill`), mirroring the dose engine design. The sequencer is selected via `pln.propSeq.sequencer` and discovered automatically. VMAT (dynamic/arc) delivery is now also implemented in this hierarchy: `matRad_PhotonLeafSequencerVMATAbstract` holds the FMO-beam gating, arc spreading and leaf-speed/dose-rate post-processing shared by any VMAT-capable sequencer, and `matRad_PhotonLeafSequencerSiochi` (the only sequencer that ever supported VMAT) derives from it. The previous `matRad_siochiLeafSequencing`/`matRad_xiaLeafSequencing`/`matRad_engelLeafSequencing` functions have been removed; `pln.propSeq.sequencer`/`pln.propOpt.runVMAT` keep working exactly as before.
- Argument order of `matRad_sequencing` changed to `matRad_sequencing(resultGUI, stf, pln, dij)`; the previous order (with `dij` and `pln` swapped) is still accepted with a deprecation warning
- The arc bookkeeping in `stf` and `apertureInfo` was renamed for readability: the `propVMAT` sub-struct became `arc`, `lim_l`/`lim_r` became `limLeft`/`limRight`, the `_I`/`_F` leaf position suffixes became `Initial`/`Final`, and `apertureInfo.constraints` became `apertureInfo.deliveryConstraints`. `matRad_upgradeApertureInfo` converts aperture info saved by an older matRad version, so stored plans keep working
- The sequencer property `sequencingLevel` was renamed to `numLevels`; the old name still works as a deprecated alias. For VMAT it is only a starting value, since the sequencer re-stratifies at increasing level counts until every fluence-optimized beam yields at least as many apertures as it has DAO control points
- `matRad_doseRecalc` was replaced by `matRad_refineApertureArc`, which separates interpolating the apertures onto a finer gantry angle grid from the forward dose calculation, and can optionally reuse an existing dij instead of recomputing one
- `matRad_calcDeliveryMetrics` moved to `planAnalysis` and is now a side-effect-free calculation returning the metrics in `resultGUI.deliveryMetrics`, instead of printing them and modifying its input
- `matRad_fluenceOptimization` no longer needs the `stf`: the fluence-optimization beam bookkeeping of a VMAT plan travels in the dij (`dij.isFMOBeam`). Passing an stf in the 4th position is still accepted with a deprecation warning
- The post-DAO prescription scaling moved from the top-level `pln.scaleDRx`/`DRx`/`RxStruct` to `pln.propOpt.scaleToPrescription`/`prescribedDose`/`prescriptionStructIx`; the old locations are honored with a deprecation warning
- The DAO and VMAT optimization problems no longer keep near-duplicate copies of the objective, gradient and Jacobian code: `matRad_OptimizationProblemVMAT` derives from `matRad_OptimizationProblemDAO` and overrides only what differs. `matRad_bixWeightAndGrad` was renamed `matRad_calcBixelWeightAndGradient` and its contract aligned with the API names
- The VMAT gantry angle helper `matRad_VMATGantryAngles` was removed; the arc angle placement it performed is part of the `PhotonVMAT` stf generator
- The Gaussian fluence smoothing that VMAT sequencing applied to every fluence map before stratifying it is now selectable via `pln.propSeq.arcFluenceSmoothing` (`'gaussian'`, the default, or `'none'`), and how surplus apertures are discarded via `pln.propSeq.apertureSelection` (`'doseAreaProduct'`, the default, or `'leastSquares'`, which keeps the apertures that best reconstruct the optimized fluence and refits their weights to it). The loop that raises the stratification level count until enough apertures exist is now bounded and reports what to change, instead of spinning forever on a fluence that is flat over its whole field
- The promotion of cst objectives/constraints to objects and the rescaling of their dose parameters to the dose per fraction, previously duplicated in `matRad_fluenceOptimization` and `matRad_directApertureOptimization`, moved into `matRad_fractionateCstFunctions`

## 3.2.3 - 2026-07-11

### Added

- Support for oxygen ions, including `oxygen_Generic` base data, the MKM biological model (`matRad_MKM` / z*-based LQ models), and oxygen handling in the particle stf generators and pencil-beam engine
- Userfolders can now also be set via environment variable `MATRAD_USERDATA`
- Documentation: Documented the userfolder feature and its usage as well as other datastructures more clearly
- Range shifter lateral scattering (`matRad_calcSigmaRashi`) is now modeled for heavier ions (helium, carbon, oxygen) in addition to protons
- Precompiled IPOPT interface binaries for Octave 8.4.0 on Linux (`ipopt.mexoct840a64`) and Windows (`ipopt.mexoct840w64`), together with an updated MinGW compilation script
- CI: added dependabot configuration for monthly, grouped GitHub Actions updates
- Octave >= 10 mex file for ipopt linked against Octave's OpenBLAS and LAPACK.

### Fixed

- possible negative doses in finesampling engine due to extrapolation in kernel interpolation
- correct parsing of all optional arguments of the `traceCube` function for `matRad_RayTracer`
- corrected an inconsistency where analytical dose calculations used an RSP cube with `ignoreOutsideDensities` applied, while MC dose engines often converted materials directly from the HU cube without the same masking. Added the variables `ignoreOutsideDensities` and `useGivenEqDensityCube` to `pln.propStf` and `pln.propDoseCalc` to handle this consistently between STF generation and dose engines. Note that `ignoreOutsideDensities` now defaults to `false` (densities outside contours are kept), whereas dose calculation and STF generation previously masked outside densities by default, so results can differ unless the option is set explicitly.
- CheckGradients option for fmincon dropped due to change in Matlab 2026
- `dij.ax` and `dij.bx` are correctly handled as cell arrays in weight initialization
- `matRad_GriddedScenariosAbstract` now correctly allows single grid point (collapse to nominal value) for a specific error type
- Range shifter lateral scattering is now correctly applied to the lateral dose kernels of the analytical pencil-beam engines (HongPB, AnalyticalPB, SubsamplingPB); previously it was computed but never added, so range shifters produced no lateral broadening (#923)
- Fixed a unit bug in `matRad_calcSigmaRashi` where the base data range (in mm) was used directly in a cm-based formula, underestimating the range shifter scattering for protons
- Correct edge case in handling of empty bixels in `matRad_ParticlePencilBeamEngineAbstract`
- CI: pinned the Octave dicom package to 0.7.2 due to a build bug in newer releases

### Changed

- New version of photons_Generic.mat basedata file can now be provided, allowing a "version" field alongside "meta" and "data" files within the machine struct. Version 2 requires correct kernel normalization (without implying a spacing in the convolution integral). photons_Generic.mat has been updated to version 2 with correct kernel normalization.
- Photon dose calculation now does not rely on hardcoded convolution resolution integral normalization of machine kernels. Assumes that old kernels use hardcoded factor of 4 for 0.5 mm resolution (1/0.5^2).
- Improved matching of RTStruct contours to ct slices in DICOM import
- CI: Octave tests now run with Octave 8.4.0 on ubuntu-24.04 (previously Octave 6.4 on ubuntu-22.04)
- CI: updated GitHub Actions to current major versions (checkout v7, upload-artifact v7, download-artifact v8, and others)
- Octave now manages mex file versions differently. Since Octave 10, mex files only link against a dedicated mex library (instead of full octave and libinterp1). The mex file checker now checks for the latest available major version build, and tries to run it.

## [3.2.2] - 2026-03-26

### Fixed

- Fixed `matRad_version` in case of tagged releases.
- Fixed documentation: small documentation correction.

## [3.2.1] - 2026-03-25

### Added

- This patch fixes a multitude of issues with the new MATLAB Desktop from R2025 and reported minor issues in import/export, helper functions, and other utility functions. Apart from that, it introduces new flexibilities under the hood, such as single precision and raytracer vectorization, that do not intentionally change outward-facing behavior and can be tested before broader behavioral changes in a future major release. The license was also changed to BSD 3-Clause.
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

## [3.2.0] - 2025-10-30

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

## [3.1.0] - 2024-11-18

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

## [2.10.1] - 2020-11-20

### Added

- Added blue/white/red difference maps to the available colormaps.
- Added the option `pln.propDoseCalc.useGivenEqDensityCube`, defaulting to `false`, to directly use the literal values from `ct.cube` and omit HU to WEQ conversion from `ct.cubeHU`.
- Added the option `pln.propDoseCalc.ignoreOutsideDensities`, defaulting to `true`, to disable or enable inclusion of WEPL outside the patient contour during ray tracing.

### Changed

- Release with small updates, cleanups, and bug fixes.
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

## [2.10.0] - 2020-06-05

### Added

- Second release of matRad. Despite major incompatibilities with "Alan", the project kept the major version number `2` to preserve a more consistent versioning scheme going forward. The team thanks all new contributing authors listed in `AUTHORS.txt`.
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

## [2.1.0] - 2016-05-23

### Added

- First official release of matRad.
- Added the IPOPT optimizer for constrained optimization.
- Added validated ray tracing.
- Added validated pencil beam particle dose calculation.
- Added validated singular value decomposed pencil beam photon dose calculation.
- Added DICOM import including dose and particle pencil beam scanning plan objects.
- Added a standalone matRad version.
- Improved the GUI workflow.

### Fixed

- Many bug fixes and many new bugs.

[Unreleased]: https://github.com/e0404/matRad/compare/v3.2.2...HEAD
[3.2.2]: https://github.com/e0404/matRad/compare/v3.2.1...v3.2.2
[3.2.1]: https://github.com/e0404/matRad/compare/v3.2.0...v3.2.1
[3.2.0]: https://github.com/e0404/matRad/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/e0404/matRad/compare/v2.10.1...v3.1.0
[2.10.1]: https://github.com/e0404/matRad/compare/v2.10.0...v2.10.1
[2.10.0]: https://github.com/e0404/matRad/compare/2.1.0...2.10.0
[2.1.0]: https://github.com/e0404/matRad/compare/initial...2.1.0
