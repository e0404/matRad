# Collapsed cone kernel generation

Scripts to generate the energy deposition kernels used by
`DoseEngines.matRad_PhotonCollapsedConeEngine` with EGSnrc/EDKnrc and to
fit the collapsed cone parameters.

The resulting data file is
`matRad/doseCalc/collapsedCone/matRad_collapsedConeKernel_6MV.mat`.
Raw EGSnrc output is not part of the repository.

## Requirements

- EGSnrc configured, with the `edknrc` user code compiled
  (`cd $EGS_HOME/edknrc && make`)
- MATLAB or Octave with matRad on the path

## Steps

1. Generate the monoenergetic input files (0.25 to 6 MeV, water,
   `521icru.pegs4dat`, source 0, 48 cones of 3.75 deg, radii from 0.05 to
   60 cm). The optional argument is the number of histories per energy
   (default 1e7):

   ```bash
   ./writeEdkInputs.sh 10000000
   ```

   The generated files are stored in `egsinp/`.

2. Run all energies in parallel (one job per energy, `nice -n 19`). Listings,
   logs and wall times (`walltime.txt`) are written to the given folder:

   ```bash
   ./runEdknrc.sh /path/outside/the/repository/edk_runs
   ```

   With 1e7 histories per energy the 9 jobs take about 2.7 h on 8 cores
   (between 0.5 h for 0.25 MeV and 2.7 h for 6 MeV).

3. Build the polyenergetic kernel and fit the collapsed cone parameters:

   ```matlab
   matRad_rc;
   addpath('tools/collapsedConeKernels');
   matRad_buildCollapsedConeKernel('/path/outside/the/repository/edk_runs', ...
       fullfile(getenv('HEN_HOUSE'), 'spectra', 'egsnrc', 'mohan6.spectrum'), [], 12, 1.5, 200);
   ```

   The last three arguments are the number of polar bins of the stored fit,
   the smallest radius included in the fit (1.5 mm, half of a 3 mm voxel) and
   the largest radius included in the fit of the scatter component (200 mm).

## Method

- `matRad_readEdkListing` reads the conical kernel table (total and primary
  energy fractions per cone and shell) from each `.egslst` file. The scatter
  component is total minus primary.
- `matRad_combineEdkSpectrum` weights the monoenergetic kernels with the TERMA
  fraction `Phi(E) * E * mu/rho(E)` of water for the Mohan 6 MV spectrum,
  interpolating linearly between the simulated energies.
- `matRad_fitCollapsedConeKernel` (in `matRad/doseCalc/collapsedCone`)
  collapses the cones into polar bins and fits, per bin,
  `dE/dr = A a exp(-a r) + B b exp(-b r)` (r in mm of water). `A` and `B` are
  the primary and scatter energy fractions of the bin from the Monte Carlo
  data; `a` and `b` are fitted to the normalized cumulative radial energy
  distribution of each component for r >= 1.5 mm. Energy deposited closer to
  the interaction point stays in the source voxel of the transport. The
  scatter is fitted only up to 200 mm: EDKnrc scores inside a 60 cm sphere,
  so the normalized cumulative scatter of the data reaches 1 at 60 cm (the
  escaping energy is missing), which a single exponential cannot follow, and
  distances of patient size matter most for the dose. The scatter build-up
  near the origin cannot be represented by positive exponentials.
- No renormalization is applied: the energy escaping the 60 cm sphere is
  reported as missing fraction, `sum(A + B) = 1 - missing`. Absolute dose
  normalization is performed by the dose engine against the machine data.
- The monoenergetic kernels are stored in the data file as well.
