#!/usr/bin/env bash
set -euo pipefail

outDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/egsinp"
histories="${1:-10000000}"
energies=(0.25 0.5 1.0 1.5 2.0 3.0 4.0 5.0 6.0)

mkdir -p "$outDir"
seed=1
for e in "${energies[@]}"; do
    tag=$(printf "%05d" "$(awk -v e="$e" 'BEGIN{printf "%d", e*1000+0.5}')")
    seed=$((seed + 1))
    cat > "$outDir/edk_water_${tag}keV.egsinp" <<INP
TITLE= Energy deposition kernel in water for ${e} MeV photons

:start I/O control:
IRESUME= first
STORE DATA ARRAYS= yes
PRINT OUT EDK FILE= no
:stop I/O control:

:start Monte Carlo inputs:
NUMBER OF HISTORIES= ${histories}
INITIAL RANDOM NO. SEEDS= ${seed}, $((seed * 97))
IFULL= ENERGY DEPOSITION KERNEL
DOPPLER BROADENING= On
:stop Monte Carlo inputs:

:start geometrical inputs:
NUMBER OF CONES= 48
ANGLES= 3.75
RADII= 0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1.0,
       1.5,2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0,20.0,
       30.0,40.0,50.0,60.0
MEDIA= H2O521ICRU;
MEDNUM= 1
START REGION= 2
STOP REGION= 1153
:stop geometrical inputs:

:start source inputs:
INCIDENT PARTICLE= photon
INCIDENT ENERGY= monoenergetic
INCIDENT KINETIC ENERGY(MEV)= ${e}
SOURCE NUMBER= 0
:stop source inputs:

:start MC transport parameter:
Global ECUT= 0.521
Global PCUT= 0.010
Global SMAX= 0.0
ESTEPE= 0.1
XImax= 0.0
Skin depth for BCA= 3
Boundary crossing algorithm= EXACT
Electron-step algorithm= PRESTA-II
Spin effects= on
Brems angular sampling= KM
Brems cross sections= BH
Bound Compton scattering= On
Pair angular sampling= Simple
Photoelectron angular sampling= On
Rayleigh scattering= On
Atomic relaxations= On
Photon cross sections= si
:stop MC transport parameter:

:start variance reduction:
ELECTRON RANGE REJECTION= on
ESAVEIN= 1.0
EXPONENTIAL TRANSFORM C= 0.0000
PHOTON FORCING= On
START FORCING= 1
STOP FORCING AFTER= 2
:stop variance reduction:

:start plot control:
PLOTTING= Off
:stop plot control:
INP
done
