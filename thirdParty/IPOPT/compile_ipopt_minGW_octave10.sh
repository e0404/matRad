# =============================================================================
# Compile the IPOPT mex interface for Octave 10 (64-bit) on Windows
# =============================================================================
#
# matRad ships the IPOPT mex interface pre-compiled for Matlab. This script
# builds it from the MSYS2/MinGW distribution that comes bundled with Octave
# for Windows, so that it can be used from within Octave.
#
# WHAT IS DIFFERENT FROM THE OLDER SCRIPTS (6.4.0 / 8.4.0 / 9.2.0):
#
#   1. IPOPT is updated from 3.12.13 to 3.14.19. This is REQUIRED, not
#      cosmetic: 64-bit integer support (--with-intsize=64) only exists in
#      IPOPT >= 3.14 and in the modern coin-or-tools ThirdParty-Mumps. The
#      3.14 series also no longer bundles the ThirdParty sources in-tree, so
#      MUMPS is built and installed separately and located via pkg-config.
#
#   2. We link against the OpenBLAS that ships WITH Octave
#      ($MINGW_PREFIX/lib/libopenblas.dll.a, resolving $MINGW_PREFIX/bin/
#      libopenblas.dll at runtime) instead of building the reference
#      COIN BLAS/LAPACK. No -lcoinblas / -lcoinlapack anymore.
#
#   3. Octave's OpenBLAS is built with 64-bit array indexing (ILP64), i.e.
#      its BLAS/LAPACK symbols expect 8-byte integer arguments. Anything that
#      calls it must therefore also pass 8-byte integers. We build BOTH MUMPS
#      and IPOPT with --with-intsize=64 (which passes -fdefault-integer-8 to
#      gfortran and defines the 64-bit integer types in the C/C++ headers).
#
#   4. METIS IS DROPPED. ThirdParty-Metis only packages METIS 4.0.3, whose
#      idx_t is a 32-bit int, and ThirdParty-Mumps refuses to build a 64-bit
#      MUMPS against a 32-bit METIS (hard #error in its configure check).
#      Building a 64-bit METIS is a separate patching exercise. MUMPS works
#      without METIS (it uses its built-in PORD / AMD orderings), so we pass
#      --without-metis. For matRad's problem sizes the ordering difference is
#      not relevant. (If you ever want METIS back, you must supply a METIS 5
#      built with IDXTYPEWIDTH=64 and point MUMPS at it -- see notes at bottom.)
#
#   5. Octave 10 links .mex files against the new liboctmex library instead of
#      liboctinterp + liboctave. `mkoctfile --mex` handles this automatically,
#      so the compile line does not need to name any Octave library.
#
# HOW TO RUN:
#   Start an Octave MinGW shell (run "cmdshell.bat" from your Octave install
#   directory, if needed as administrator) and execute this script from there.
#   The script stops at the first error so failures are easy to locate; just
#   re-run after fixing. Intermediate build trees live under the current
#   working directory (./ipopt-build).
# =============================================================================

set -e

# ---------------------------------------------------------------------------
# 0) Tools and configuration
# ---------------------------------------------------------------------------
pacman -Sy
pacman -S --noconfirm --needed which

echo "check_certificate = off" >> ~/.wgetrc

# Pinned versions (bump here if you want to iterate)
IPOPT_VERSION=3.14.19
MUMPS_BRANCH=releases/3.0.12          # coin-or-tools/ThirdParty-Mumps (MUMPS 5.x)
MEXIPOPT_BRANCH=1.1.6                 # ebertolazzi/mexIPOPT
METIS_VERSION=5.1.0                   # KarypisLab METIS, built with 64-bit idx_t

# Build MUMPS with OpenMP (multithreaded factorization)?  Default 0 (sequential).
# Sequential is recommended for a DISTRIBUTABLE mex: a threaded MUMPS nests with
# Octave's threaded OpenBLAS and oversubscribes cores -- you get the runtime
# "OpenBLAS warning: precompiled NUM_THREADS exceeded" and a per-iteration
# slowdown, avoidable only by tuning OMP_NUM_THREADS at runtime. A sequential
# MUMPS sidesteps all of that (OpenBLAS still threads the dense algebra). Set to
# 1 only if you specifically want threaded MUMPS and will cap threads yourself.
# NOTE: either way the build needs omp_lib.mod present, because MUMPS' sources
# do an unconditional `USE OMP_LIB`; Stage 0b restores it regardless.
MUMPS_OPENMP=0

# Octave's bundled ILP64 LAPACK + OpenBLAS. $MINGW_PREFIX is /mingw64 in the
# Octave MinGW shell. NOTE: Octave's libopenblas is BLAS-only here (it does NOT
# contain LAPACK); LAPACK ships as a separate liblapack.dll.a in
# $MINGW_PREFIX/lib. Both are built with 64-bit integer indexing (ILP64) --
# that is what Octave itself uses -- so they match our --with-intsize=64 build.
# -llapack must come before -lopenblas (LAPACK depends on BLAS).
export LAPACK_LFLAGS="-L$MINGW_PREFIX/lib -llapack -lopenblas"

# Force libtool to use the clean MSYS linker path. Left to auto-detection,
# libtool asks the compiler for `ld` and gets the Windows-style path
# "C:/Program Files/GNU Octave/.../ld.exe"; the spaces in it get mangled into
# newlines and end up splitting the generated Makefile (=> "missing separator"
# errors). $MINGW_PREFIX/bin/ld is the same linker with a space-free path.
export lt_cv_path_LD="$MINGW_PREFIX/bin/ld"

# Working / install directories
WORKDIR="$(pwd)/ipopt-build"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

mkdir -p install
export PREFIX="$(cd install && pwd)"
# ThirdParty-Mumps installs its coinmumps.pc / ipopt.pc here so configure and
# the final mex link can find everything via pkg-config.
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"

# Sanitize generated .pc files. autoconf detects gfortran's runtime library
# directory, which on this Octave install is a Windows-style path with spaces
# ("C:/Program Files/GNU Octave/.../lib"); it gets truncated at the first space
# into a broken "-LC:/Program" token that pollutes coinmumps.pc / ipopt.pc.
# Everything those -L paths provide (libgfortran/libquadmath/libgomp/...) is
# also reachable through -L$MINGW_PREFIX/lib, which is already on the link line,
# so we simply drop any drive-letter -L path from the .pc files.
sanitize_pc() {
  for f in "$PREFIX"/lib/pkgconfig/*.pc ; do
    [ -f "$f" ] && sed -i -E 's, -L[A-Za-z]:[^ ]*,,g' "$f"
  done
}

echo "==================================================================="
echo " Build prefix : $PREFIX"
echo " LAPACK/BLAS  : $LAPACK_LFLAGS"
echo " IPOPT        : $IPOPT_VERSION"
echo " MUMPS branch : $MUMPS_BRANCH (--with-intsize=64, METIS $METIS_VERSION 64-bit, OpenMP=$MUMPS_OPENMP)"
echo "==================================================================="

# ---------------------------------------------------------------------------
# 0b) Restore the gfortran OpenMP module (omp_lib.mod)
# ---------------------------------------------------------------------------
# Octave's Windows build (MXE) ships gfortran WITHOUT the Fortran OpenMP
# support files: its .../finclude directory contains only the ieee_* modules,
# no omp_lib.mod / omp_lib_kinds.mod / omp_lib.h (the C side, omp.h + libgomp,
# is present). Current MUMPS (>= 5.8.x) has an UNCONDITIONAL `USE OMP_LIB`, so
# the build cannot proceed without that module regardless of --enable-openmp.
# We restore it from the *same GCC version's* MSYS2 gcc-fortran package (.mod
# files are tied to the GCC major version, so the version must match exactly;
# we detect it rather than hardcoding). This writes into Octave's install tree,
# so run the script from an Administrator MinGW shell.
echo "==== STAGE 0b: ensure gfortran OpenMP module (omp_lib.mod) ===="
FINCLUDE="$(gfortran -print-file-name=finclude)"
if [ -f "$FINCLUDE/omp_lib.mod" ]; then
  echo "omp_lib.mod already present in $FINCLUDE -- skipping restore."
else
  GFVER="$(gfortran -dumpfullversion)"     # e.g. 14.2.0
  GFTRIPLE="$(gcc -dumpmachine)"           # e.g. x86_64-w64-mingw32
  MSYS_REPO="https://repo.msys2.org/mingw/mingw64"
  echo "omp_lib.mod missing; restoring from MSYS2 gcc-fortran $GFVER ($GFTRIPLE)."
  # Pick the newest package revision matching this exact gfortran version.
  PKG=$(wget --no-check-certificate -qO- "$MSYS_REPO/" \
        | grep -oE "mingw-w64-x86_64-gcc-fortran-${GFVER}-[0-9]+-any\.pkg\.tar\.zst" \
        | sort -V | tail -n1)
  if [ -z "$PKG" ]; then
    echo "ERROR: no MSYS2 gcc-fortran package found for gfortran $GFVER at $MSYS_REPO." >&2
    echo "       Fetch omp_lib.mod / omp_lib_kinds.mod for GCC $GFVER manually into:" >&2
    echo "       $FINCLUDE" >&2
    exit 1
  fi
  OMPD="$WORKDIR/ompmod"
  rm -rf "$OMPD" && mkdir -p "$OMPD" && cd "$OMPD"
  wget --no-check-certificate "$MSYS_REPO/$PKG"
  # Decompress explicitly with zstd (GNU tar needs the zstd CLI for .zst, which
  # is not guaranteed; pacman itself depends on zstd so this is cheap/safe).
  pacman -S --noconfirm --needed zstd
  zstd -dqf "$PKG"                          # -> ${PKG%.zst} (the .tar)
  # Extract the whole finclude dir (grabs omp_lib.mod, omp_lib_kinds.mod, .h).
  tar -xf "${PKG%.zst}" "mingw64/lib/gcc/$GFTRIPLE/$GFVER/finclude"
  cp "$OMPD/mingw64/lib/gcc/$GFTRIPLE/$GFVER/finclude/"omp_lib* "$FINCLUDE"/
  cd "$WORKDIR"
  # Verify a minimal `use omp_lib` now compiles.
  printf '      program t\n      use omp_lib\n      print *, omp_get_max_threads()\n      end\n' > "$OMPD/omptest.f"
  if gfortran -fopenmp "$OMPD/omptest.f" -o "$OMPD/omptest.exe" ; then
    echo "omp_lib.mod restored into $FINCLUDE and verified."
  else
    echo "ERROR: omp_lib.mod restore failed -- 'use omp_lib' still does not compile." >&2
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# 0c) Build METIS 5.1.0 with 64-bit integer indices (idx_t = int64)
# ---------------------------------------------------------------------------
# A 64-bit MUMPS can only use METIS if METIS' idx_t is also 64-bit: MUMPS'
# configure has a hard #error unless metis.h reports IDXTYPEWIDTH == 64. The
# prebuilt MSYS2 metis package is 32-bit, and ThirdParty-Metis only packages
# the ancient METIS 4, so we build METIS 5.1.0 ourselves with IDXTYPEWIDTH=64.
# METIS 5 uses CMake + a bundled GKlib that needs a few MinGW fixes (the same
# ones MSYS2 applies); we patch them inline and build only the library target
# (the gpmetis/ndmetis command-line programs are not needed and are the most
# troublesome to build on MinGW). Installed into $PREFIX alongside MUMPS/IPOPT.
echo "==== STAGE 0c: METIS $METIS_VERSION (64-bit idx_t) ===="
# CMake is required to build METIS. Octave's Windows bundle already ships cmake
# in $MINGW_PREFIX/bin, so do NOT pacman-install it: Octave's /mingw64 tree is
# not pacman-managed, so pulling cmake's dependency chain both conflicts with
# the existing (untracked) files and risks upgrading the gcc runtime away from
# the version we matched omp_lib.mod to. Use Ninja if it is on PATH, otherwise
# fall back to MSYS Makefiles (the MinGW shell always provides `make`).
if ! command -v cmake >/dev/null 2>&1 ; then
  echo "ERROR: cmake not found on PATH (Octave normally ships it in $MINGW_PREFIX/bin)." >&2
  echo "       Install a standalone CMake and re-run, or drop METIS (see note at end)." >&2
  exit 1
fi
if command -v ninja >/dev/null 2>&1 ; then METIS_GEN="Ninja"; else METIS_GEN="MSYS Makefiles"; fi
echo "Using CMake $(cmake --version | head -n1) with generator: $METIS_GEN"
# Download and patch the source only on a fresh extract, so the (non-idempotent)
# sed patches are never re-applied on a re-run of the script.
if [ ! -d "metis-$METIS_VERSION" ]; then
  wget --no-check-certificate \
    "https://papers.karypis.org/glaros/files/sw/metis/metis-$METIS_VERSION.tar.gz"
  tar -xzf "metis-$METIS_VERSION.tar.gz"
  cd "metis-$METIS_VERSION"
  # --- MinGW patches (equivalent to MSYS2's mingw-w64-metis patches) ---
  # 1) GKlib includes <sys/resource.h>, which MinGW lacks.
  sed -i 's|\(#include <sys/resource.h>\)|#ifndef __MINGW32__\n\1\n#endif|' GKlib/gk_arch.h
  # 2) GKlib's getopt prototypes use reserved names __argc/__argv (real symbols
  #    on MinGW); rename them so the headers compile.
  sed -i 's/__argc/gk_argc/g; s/__argv/gk_argv/g' GKlib/gk_getopt.h
  # 3) Make GKLIB_PATH absolute so out-of-tree CMake builds find GKlib.
  sed -i 's|set(GKLIB_PATH "GKlib"|set(GKLIB_PATH "${CMAKE_SOURCE_DIR}/GKlib"|' CMakeLists.txt
  # --- select 64-bit integer indices ---
  sed -i -E 's/^#define IDXTYPEWIDTH[[:space:]]+32/#define IDXTYPEWIDTH 64/' include/metis.h
  cd "$WORKDIR"
fi
cd "metis-$METIS_VERSION"
grep -E '^#define IDXTYPEWIDTH' include/metis.h   # sanity: must show 64

# Configure + build ONLY the library (CMAKE_POLICY_VERSION_MINIMUM lets modern
# CMake accept METIS' ancient cmake_minimum_required(VERSION 2.8)).
rm -rf build
cmake -S . -B build -G "$METIS_GEN" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Release \
  -DSHARED=OFF \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build --target metis
# Install just the pieces MUMPS needs (libmetis.a already contains GKlib).
# METIS is built before MUMPS, so $PREFIX/lib and /include may not exist yet.
mkdir -p "$PREFIX/lib" "$PREFIX/include"
cp "$(find build -name 'libmetis.a' | head -1)" "$PREFIX/lib/"
cp include/metis.h "$PREFIX/include/"
cd "$WORKDIR"
# METIS link/compile flags for the MUMPS configure below.
# -Wno-error=incompatible-pointer-types: ThirdParty-Mumps' metis availability
# check calls METIS_NodeND((int*)0, ...) with int* regardless of idx_t width.
# With 64-bit idx_t (long long*) that is an incompatible pointer type, which
# GCC 14 makes a hard error by DEFAULT, so the otherwise-harmless check fails
# to compile. Downgrading it back to a warning lets the check (and the build)
# succeed; MUMPS' real metis interface uses the correct idx_t types.
export METIS_CFLAGS="-I$PREFIX/include -Wno-error=incompatible-pointer-types"
export METIS_LFLAGS="-L$PREFIX/lib -lmetis"

# ---------------------------------------------------------------------------
# 1) Build MUMPS (64-bit integers; Octave OpenBLAS/LAPACK; 64-bit METIS)
# ---------------------------------------------------------------------------
echo "==== STAGE 1: ThirdParty-Mumps (64-bit) ===="
if [ ! -d ThirdParty-Mumps ]; then
  git clone --depth 1 --branch "$MUMPS_BRANCH" \
    https://github.com/coin-or-tools/ThirdParty-Mumps.git
fi
cd ThirdParty-Mumps

# Fetch the MUMPS source (the default coin-or-tools.github.io mirror is used).
./get.Mumps

# Insurance against the classic "cannot find -lquadmath" from gfortran: make a
# copy of libquadmath where a bare `gcc` search path expects it. Harmless if
# unused by the autotools build.
export GCC_VERSION=`gcc --version | head -n1 | cut -d" " -f3`
mkdir -p /usr/lib/gcc/$MINGW_CHOST/$GCC_VERSION/
cp $MINGW_PREFIX/lib/gcc/$MINGW_CHOST/$GCC_VERSION/libquadmath.* \
   /usr/lib/gcc/$MINGW_CHOST/$GCC_VERSION/ 2>/dev/null || true

rm -rf build && mkdir build && cd build
# OpenMP toggle (see MUMPS_OPENMP above). Without -fopenmp the `!$omp` regions
# become comments (sequential MUMPS); the unconditional `USE OMP_LIB` still
# compiles because Stage 0b put omp_lib.mod on the module search path.
if [ "$MUMPS_OPENMP" = "1" ]; then
  MUMPS_OPENMP_CONFIG="--enable-openmp"
  echo "Building MUMPS WITH OpenMP (threaded)."
else
  MUMPS_OPENMP_CONFIG="--disable-openmp"
  echo "Building MUMPS WITHOUT OpenMP (sequential)."
fi
../configure --prefix="$PREFIX" \
  --disable-shared --enable-static \
  --disable-dependency-tracking \
  $MUMPS_OPENMP_CONFIG \
  --with-intsize=64 \
  --with-metis \
  --with-metis-cflags="$METIS_CFLAGS" \
  --with-metis-lflags="$METIS_LFLAGS" \
  --with-lapack-lflags="$LAPACK_LFLAGS" \
  ADD_FCFLAGS=-fallow-argument-mismatch
make
make install
sanitize_pc   # strip broken drive-letter -L paths from coinmumps.pc
cd "$WORKDIR"

# ---------------------------------------------------------------------------
# 2) Build IPOPT (64-bit integers, OpenBLAS, MUMPS via pkg-config)
# ---------------------------------------------------------------------------
echo "==== STAGE 2: IPOPT $IPOPT_VERSION (64-bit) ===="
if [ ! -d Ipopt ]; then
  git clone --depth 1 --branch "releases/$IPOPT_VERSION" \
    https://github.com/coin-or/Ipopt.git
fi
cd Ipopt
rm -rf build && mkdir build && cd build
# --with-mumps is on by default and located through PKG_CONFIG_PATH.
../configure --prefix="$PREFIX" \
  --disable-shared --enable-static \
  --disable-dependency-tracking \
  --with-intsize=64 \
  --with-lapack-lflags="$LAPACK_LFLAGS" \
  --disable-java
make
make install
sanitize_pc   # strip broken drive-letter -L paths from ipopt.pc / coinmumps.pc
cd "$WORKDIR"

# Sanity check: pkg-config must now describe a usable IPOPT.
echo "---- pkg-config ipopt cflags: ----"; pkg-config --cflags ipopt
echo "---- pkg-config ipopt static libs: ----"; pkg-config --static --libs ipopt

# ---------------------------------------------------------------------------
# 3) Build the mex interface
# ---------------------------------------------------------------------------
echo "==== STAGE 3: mexIPOPT ===="
if [ ! -d mexIPOPT ]; then
  git clone --depth 1 --branch "$MEXIPOPT_BRANCH" \
    https://github.com/ebertolazzi/mexIPOPT.git
fi

# `mkoctfile --mex` links against liboctmex automatically on Octave 10.
# The include dir and the full static dependency chain (IPOPT + coinmumps +
# OpenBLAS + pord ...) come straight from pkg-config so we don't hardcode
# paths or the include/coin-or vs include/coin naming.
# -lgomp: MUMPS references the OpenMP runtime symbols (via `USE OMP_LIB`) even
# in the sequential build, so link libgomp explicitly. Harmless in the OpenMP
# build too, where pkg-config already pulls it in.
mkoctfile --mex \
  -ImexIPOPT/toolbox/src \
  $(pkg-config --cflags ipopt) \
  mexIPOPT/toolbox/src/ipopt.cc \
  mexIPOPT/toolbox/src/IpoptInterfaceCommon.cc \
  -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N \
  $(pkg-config --static --libs ipopt) \
  -lgfortran -lquadmath -lgomp

echo "==================================================================="
echo " Done. If the build succeeded, ipopt.mex is in $WORKDIR."
echo " matRad ships the archive copy under a MAJOR-version name so any Octave"
echo " >= 10 picks the newest compatible one (see matRad_checkMexFileExists):"
echo "   copy it next to that function as  ipopt.mexoct10w64  (Octave 10.x, win64)."
echo "==================================================================="

# ---------------------------------------------------------------------------
# Cleanup (optional)
# ---------------------------------------------------------------------------
echo "Do you wish to remove the IPOPT compilation directories (ThirdParty-Mumps, Ipopt)?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) rm -rf "$WORKDIR/ThirdParty-Mumps" "$WORKDIR/Ipopt"; break;;
        No ) exit;;
    esac
done

# =============================================================================
# NOTE on METIS: METIS 5.1.0 is built above (Stage 0c) with IDXTYPEWIDTH=64 so
# that it matches the 64-bit MUMPS. MUMPS' configure verifies at compile time
# that metis.h reports IDXTYPEWIDTH == 64 and aborts otherwise. If you ever want
# to drop METIS again (MUMPS then falls back to its built-in PORD/AMD ordering),
# skip Stage 0c and replace the three --with-metis* options in the MUMPS
# configure with a single --without-metis. -lmetis reaches the final mex link
# automatically via coinmumps.pc (pkg-config --static), so no manual edit of the
# mkoctfile line is needed either way.
# =============================================================================
