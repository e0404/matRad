# The following has been tested in Ubuntu 22.04 with gcc 
# Note that this uses a reference BLAS and LAPACK and links everything statically. 
# For better optimizer performence, consider using optimized implementations (e.g. OpenBLAS) and an updated IPOPT and use this script as a guideline.
# You could also link against openblas and lapack as octave uses those libraries anyways
echo "check_certificate = off" >> ~/.wgetrc

# Run the following commands to create directories and get the IPOPT source. 
mkdir ipopt && cd ipopt
export IPOPTINSTALLDIR=`pwd`
cd ..

mkdir IpoptMUMPS && cd IpoptMUMPS
export IPOPTDIR=`pwd`

wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.13.tgz
tar -zxvf Ipopt-3.12.13.tgz
mv Ipopt-3.12.13/* ./

cd $IPOPTDIR/ThirdParty/Blas
# First we need to replace the url as the version can no longer be downloaded
sed -i 's,http://www.netlib.org/blas/,http://coin-or-tools.github.io/ThirdParty-Blas/,g' get.Blas
./get.Blas
cd $IPOPTDIR/ThirdParty/Lapack
./get.Lapack
cd $IPOPTDIR/ThirdParty/Metis
sed -i 's,http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/,http://coin-or-tools.github.io/ThirdParty-Metis/,g' get.Metis
./get.Metis
cd $IPOPTDIR/ThirdParty/Mumps
# First we need to replace the url as the version can no longer be downloaded
sed -i 's,http://mumps.enseeiht.fr/,http://coin-or-tools.github.io/ThirdParty-Mumps/,g' get.Mumps
./get.Mumps

cd $IPOPTDIR


# Now lets start the build process
mkdir build
cd build

../configure --prefix=$IPOPTINSTALLDIR --disable-shared --enable-static ADD_FFLAGS="-fallow-argument-mismatch -fPIC" ADD_FCLAGS="-fallow-argument-mismatch -fPIC" ADD_CXXFLAGS="-fPIC" ADD_CFLAGS="-fPIC"
make
make install

cd ../..

# If everything worked, you should see some (static) libraries when doing ls /usr/local/lib 
# we can get the mex interface

git clone --depth 1 --branch 1.1.4 https://github.com/ebertolazzi/mexIPOPT

# and compile it.
mkoctfile --mex -ImexIPOPT/toolbox/src -I$IPOPTINSTALLDIR/include/coin mexIPOPT/toolbox/src/ipopt.cc mexIPOPT/toolbox/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -lipopt -lcoinmumps -lcoinmetis -lcoinlapack -lcoinblas -lgfortran -L$IPOPTINSTALLDIR/lib

echo "Do you wish to remove the IPOPT compilation directories?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) rm -rf $IPOPTDIR; rm -rf mexIPOPT; break;;
        No ) exit;;
    esac
done


