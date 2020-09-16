# How to compile the IPOPT mex interface for Octave 5.2.0 (64-bit) in Windows
# matRad only includes the IPOPT mex interface compiled for Matlab. It is also possible to compile the interface from the MSYS/MinGW distribution included in Octave for Windows.
# The following has been tested for Octave 5.2.0 in 64 bit version to allow 64-bit algebra. Start this script from a octave mingw shell. You can open such a shell by running "cmdshell.bat" from your Octave install directory

pacman -Sy 
pacman -S --noconfirm --needed wget which git

# Run the following commands to create directories and get the IPOPT source. We dont need lapack and blas, since Octave comes with lapack and openblas

mkdir ipopt && cd ipopt
export IPOPTINSTALLDIR=`pwd`
cd ..

mkdir IpoptMUMPS && cd IpoptMUMPS
export IPOPTDIR=`pwd`

wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.13.tgz
tar -zxvf Ipopt-3.12.13.tgz
mv Ipopt-3.12.13/* ./

cd $IPOPTDIR/ThirdParty/Blas
./get.Blas
cd $IPOPTDIR/ThirdParty/Lapack
./get.Lapack
cd $IPOPTDIR/ThirdParty/Metis
./get.Metis
cd $IPOPTDIR/ThirdParty/Mumps
./get.Mumps

cd $IPOPTDIR
#- Now we need a workaround for a problem regarding the quadmath library which cannot be found when compiling MUMPS. Run the following command
export GCC_VERSION=`gcc --version | head -n1 | cut -d" " -f3`
mkdir -p /usr/lib/gcc/$MINGW_CHOST/$GCC_VERSION/
cp $MINGW_PREFIX/lib/gcc/$MINGW_CHOST/$GCC_VERSION/libquadmath.* /usr/lib/gcc/$MINGW_CHOST/$GCC_VERSION/


# Now lets start the build process

mkdir build
cd build

../configure --prefix=$IPOPTINSTALLDIR --disable-shared --enable-static 
make
make install

cd ../..

# If everything worked, you should see some (static) libraries when doing ls /usr/local/lib 
# we can get the mex interface

git clone https://github.com/ebertolazzi/mexIPOPT

# and compile it.
mkoctfile --mex -ImexIPOPT/src -I$IPOPTINSTALLDIR/include/coin mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -lipopt -lcoinmumps -lcoinmetis -lcoinlapack -lcoinblas -lgfortran -L$IPOPTINSTALLDIR/lib

echo "Do you wish to remove the IPOPT compilation directories?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) rm -rf $IPOPTDIR; rm -rf mexIPOPT; break;;
        No ) exit;;
    esac
done


