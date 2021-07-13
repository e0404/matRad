#!/bin/bash
# use as follows:
# put install script in folder together with Geant4Headers*.zip and the topas tar package
# run "sudo ./installTopas.sh"

install_path="$( cd "$( dirname "$0" )" && pwd )"

# Create folder structure
mkdir $install_path/topas $install_path/topas_extensions

# Install prerequisites
apt install -y libexpat1-dev
apt install -y libgl1-mesa-dev
apt install -y libglu1-mesa-dev
apt install -y libxt-dev
apt install -y xorg-dev
apt install -y build-essential

# Install TOPAS
tar -zxvf $install_path/topas_3_*.tar.gz -C $install_path
unzip $install_path/Geant4Header*.zip -d $install_path/topas

# Install extensions
wget -O $install_path/topas_extensions/extensions.zip https://github.com/topasmc/extensions/archive/master.zip
unzip $install_path/topas_extensions/extensions.zip -d $install_path/topas_extensions
mv $install_path/topas_extensions/extensions-master/* $install_path/topas_extensions/
rm -r $install_path/topas_extensions/extensions-master
rm $install_path/topas_extensions/extensions.zip

# Install Geant4 Data
#./installGeant4 install_path

# Build
cd $install_path/topas
cmake -DTOPAS_EXTENSIONS_DIR=$install_path/topas_extensions
make

echo "Installer finished."
