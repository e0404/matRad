#!/usr/bin/env bash

sudo chmod +x runtests.sh
sudo chmod +x submodules/MCsquare/MCsquare_linux

sudo apt-get update -qq
sudo apt-get install gdb # to capture backtrace of eventual failures
# - sudo apt-get purge libopenblas-base # fixes PPA Octave 4.0 crash on Travis

# add PPA repo to use latest GNU Octave release
sudo add-apt-repository ppa:octave/stable --yes
sudo apt-get update --yes --force-yes

# dependencies
sudo apt-get install octave --yes --force-yes
sudo apt-get install liboctave-dev --yes --force-yes
sudo apt-get install git --yes --force-yes
export MATRAD=`pwd`

# IpOpt in Ubuntu repo is too old
#sudo apt-get install coinor-libipopt-dev --yes --force-yes
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

mkdir $IPOPTDIR/build
cd $IPOPTDIR/build

$IPOPTDIR/configure --prefix=/usr/local
sudo make
sudo make test
sudo make install
sudo ldconfig


# clone repository with IpOpt Interface
git clone https://github.com/ebertolazzi/mexIPOPT

# compile interface
mkoctfile --mex -ImexIPOPT/src -I/usr/local/include/coin mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -Wl,--no-undefined -lipopt -lcoinmumps -lcoinlapack -lcoinblas -L/usr/local/lib -std=gnu++11
sudo cp ipopt* $MATRAD/optimization
cd $MATRAD