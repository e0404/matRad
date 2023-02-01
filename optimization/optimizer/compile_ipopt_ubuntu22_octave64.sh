mkdir IpoptMUMPS
export IPOPTDIR=`pwd`

mkdir install
cd install
export IPOPTINSTALLDIR=`pwd`
cd ..

mkdir source
cd source
export IPOPTSRCDIR=`pwd`

wget --no-check-certificate https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.13.tgz
tar -zxvf Ipopt-3.12.13.tgz
mv Ipopt-3.12.13/* ./

# Here we can also download the new Ipopt and link against lapack / blas?
# git clone https://github.com/coin-or/Ipopt
# git clone https://github.com/coin-or-tools/ThirdParty-Mumps
# For both we do
# configure --enable-static --disable-shared ADD_CXXFLAGS="-fPIC" ADD_CFLAGS="-fPIC" ADD_FFLAGS="-fallow-argument-mismatch -fPIC" ADD_FCFLAGS="-fallow-argument-mismatch -fPIC"

cd $IPOPTSRCDIR/ThirdParty/Blas
./get.Blas
cd $IPOPTSRCDIR/ThirdParty/Lapack
./get.Lapack
cd $IPOPTSRCDIR/ThirdParty/Metis
./get.Metis
cd $IPOPTSRCDIR/ThirdParty/Mumps
# First we need to replace the url as the version can no longer be downloaded
sed -i 's,http://mumps.enseeiht.fr/,http://coin-or-tools.github.io/ThirdParty-Mumps/,g' get.Mumps
./get.Mumps

cd $IPOPTDIR

# Now lets start the build process
mkdir build
cd build
export IPOPTBUILDDIR=`pwd`

$IPOPTSRCDIR/configure --prefix=$IPOPTINSTALLDIR --disable-shared --enable-static

make
make install


