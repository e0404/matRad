#!/usr/bin/env bash

sudo chmod +x runtests.sh
sudo chmod +x submodules/MCsquare/MCsquare_linux


# dependencies
sudo zypper install octave
sudo zypper install octave-devel 
sudo zypper install git
sudo zypper install patch
export MATRAD=`pwd`

# IpOpt in Ubuntu repo is too old
# Compiling IpOpt
cd ..
mkdir IpoptMUMPS
cd IpoptMUMPS
wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.10.tgz
tar -zxvf Ipopt-3.12.10.tgz
cd Ipopt-3.12.10/
export IPOPTDIR=`pwd`

cd $IPOPTDIR/ThirdParty/Blas
./get.Blas
cd $IPOPTDIR/ThirdParty/Lapack
./get.Lapack
cd $IPOPTDIR/ThirdParty/Metis
./get.Metis
cd $IPOPTDIR/ThirdParty/Mumps
./get.Mumps

# Apply a patch to get around a bug in MUMPS
patch -p0 < mumps.patch
patch -p0 < mumps_mpi.patch
mv MUMPS/libseq/mpi.h MUMPS/libseq/mumps_mpi.h

mkdir $IPOPTDIR/build
cd $IPOPTDIR/build

$IPOPTDIR/configure --prefix=/usr/local
sudo make
sudo make test
sudo make install
sudo ldconfig

sudo zypper install gnuplot
# clone repository with IpOpt Interface
git clone https://github.com/ebertolazzi/mexIPOPT

# compile interface
mkoctfile --mex -ImexIPOPT/src -I/usr/local/include/coin mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -lipopt -lcoinmumps -lcoinlapack -lcoinblas -L/usr/local/lib -std=gnu++11
sudo cp ipopt* $MATRAD/optimization
cd $MATRAD