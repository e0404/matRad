#!/usr/bin/env bash

sudo chmod +x runtests.sh
sudo chmod +x MCsquare/bin/MCsquare_linux

sudo apt-get update -qq --yes --force-yes
sudo apt-get install gdb --yes --force-yes # to capture backtrace of eventual failures
# - sudo apt-get purge libopenblas-base --yes --force-yes # fixes PPA Octave 4.0 crash on Travis

# add PPA repo to use latest GNU Octave release
sudo add-apt-repository ppa:octave/stable --yes
sudo apt-get update --yes --force-yes

# dependencies
sudo apt-get install octave --yes --force-yes
sudo apt-get install liboctave-dev --yes --force-yes
sudo apt-get install libgdcm2-dev --yes --force-yes #for the octave dicom package
octave --no-gui --eval "pkg install -forge dicom"
octave --no-gui --eval "pkg install -forge nan"
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

sudo apt-get install gnuplot-x11 --yes --force-yes
# clone repository with IpOpt Interface
git clone https://github.com/ebertolazzi/mexIPOPT
cd mexIPOPT
git checkout 2d72b5d1b75ece34fb78d12c54116e33efdc6d42
cd ..

# compile interface
mkoctfile --mex -ImexIPOPT/toolbox/src -I/usr/local/include/coin mexIPOPT/toolbox/src/ipopt.cc mexIPOPT/toolbox/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -Wl,--no-undefined -lipopt -lcoinmumps -lcoinlapack -lcoinblas -L/usr/local/lib -std=gnu++11
sudo cp ipopt* $MATRAD/optimization
cd $MATRAD