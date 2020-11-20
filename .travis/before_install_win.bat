@echo off
set MATRAD=%~dp0

REM Compiling Ipopt
cd ..
md Ipopt
cd Ipopt

git clone https://github.com/wahln/IPOPT.git --branch eigenInterface source
md build
md install
cd build

cmake -S ../source/ -D INSTALL_PREFIX=../install/

make
make install

cd ..
cd ..


REM clone repository with IpOpt Interface
git clone https://github.com/ebertolazzi/mexIPOPT

REM compile interface
mkoctfile --mex -ImexIPOPT/src -IIpopt/install/include/ mexIPOPT/src/ipopt.cc mexIPOPT/src/IpoptInterfaceCommon.cc -v -DMATLAB_MEXFILE -DHAVE_CSTDDEF -DIPOPT_INTERFACE_MISSING_COPY_N -Wl,-undefined,error -lipopt -lblas -lf2c -llapack -LIpopt/install/lib -std=gnu++11
xcopy /s ipopt* %MATRAD%\optimization
cd $MATRAD