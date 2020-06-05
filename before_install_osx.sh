#!/usr/bin/env bash

sudo chmod +x runtests.sh
sudo chmod +x submodules/MCsquare/MCsquare_mac

export MATRAD=`pwd`


# Compiling IpOpt
cd ..
mkdir IpoptMUMPS
cd IpoptMUMPS
curl -O https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.10.tgz
tar -zxvf Ipopt-3.12.10.tgz
cd Ipopt-3.12.10/
export IPOPTDIR=`pwd`

export LIBTOOL=`which glibtool`
export LIBTOOLIZE=`which glibtoolize`

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
make &> makeIpopt.log
sudo make test
sudo make install
sudo ldconfig


# clone repository with IpOpt Interface
git clone https://github.com/ebertolazzi/mexIPOPT

# compile interfaces
#-Wl,-$MKLROOT/lib/intel64/libmkl_intel_ilp64.a -$MKLROOT/lib/intel64/libmkl_gnu_thread.a -$MKLROOT/lib/intel64/libmkl_core.a -lgomp -lpthread -lm -ldl
mkoctfile --mex -ImexIPOPT/src -I/usr/local/include/coin mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -lipopt -lcoinmumps -lcoinmetis -lcoinlapack -lcoinblas  -lgfortran -L/usr/local/lib
mkoctfile --mex -ImexIPOPT/src -I/usr/local/include/coin mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -Wl,-undefined,error -lipopt -lcoinmumps -lcoinlapack -lcoinblas -L/usr/local/lib -std=gnu++11
sudo cp ipopt* $MATRAD/optimization
cd $MATRAD