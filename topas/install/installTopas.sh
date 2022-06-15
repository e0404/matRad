#!/bin/bash
# Install-script for TOPAS 3.7
# use as follows:
# put install script in folder together with the topas tar package (only home directory works apparently)
# run "sudo ./installTopas.sh"

install_path="$( cd "$( dirname "$0" )" && pwd )"
cat topas_*_debian9.tar.gz.part_* > topas_debian9.tar.gz

# Create folder structure
mkdir $install_path/topas $install_path/topas_extensions

# Install prerequisites
apt install -y libexpat1-dev
apt install -y libgl1-mesa-dev
apt install -y libglu1-mesa-dev
apt install -y libxt-dev
apt install -y xorg-dev
apt install -y build-essential
apt install -y libharfbuzz-dev

# Install TOPAS
apt install -y unzip
tar -zxvf $install_path/topas_debian9.tar.gz -C $install_path

# Install extensions
wget -O $install_path/topas_extensions/extensions.zip https://github.com/topasmc/extensions/archive/master.zip
unzip $install_path/topas_extensions/extensions.zip -d $install_path/topas_extensions/
mv $install_path/topas_extensions/extensions-master/* $install_path/topas_extensions/
rm -r $install_path/topas_extensions/extensions-master
rm $install_path/topas_extensions/extensions.zip

# Install Geant4 Data
unzip $install_path/topas/Geant4Header*.zip -d $install_path/topas
mkdir $install_path/G4Data
cd $install_path/G4Data
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4NDL.4.6.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4EMLOW.7.9.1.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4PhotonEvaporation.5.5.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4RadioactiveDecay.5.4.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4PARTICLEXS.2.1.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4SAIDDATA.2.0.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4ABLA.3.1.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4INCL.1.0.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4PII.1.3.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4ENSDFSTATE.2.2.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4RealSurface.2.1.1.tar.gz
wget -nc -4 http://geant4-data.web.cern.ch/geant4-data/datasets/G4TENDL.1.3.2.tar.gz

tar -zxf G4NDL.4.6.tar.gz
tar -zxf G4EMLOW.7.9.1.tar.gz
tar -zxf G4PhotonEvaporation.5.5.tar.gz
tar -zxf G4RadioactiveDecay.5.4.tar.gz
tar -zxf G4PARTICLEXS.2.1.tar.gz
tar -zxf G4SAIDDATA.2.0.tar.gz
tar -zxf G4ABLA.3.1.tar.gz
tar -xzf G4INCL.1.0.tar.gz
tar -zxf G4PII.1.3.tar.gz
tar -zxf G4ENSDFSTATE.2.2.tar.gz
tar -zxf G4RealSurface.2.1.1.tar.gz
tar -zxf G4TENDL.1.3.2.tar.gz

export TOPAS_G4_GATA_DIR= $install_path/G4Data

# Build
cd $install_path/topas
apt install -y cmake
cmake -DTOPAS_EXTENSIONS_DIR=$install_path/topas_extensions
make

echo "Installer finished."
