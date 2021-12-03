#!/bin/bash
mkdir $1/G4Data
unzip $1/Geant4Header*.zip -d $1/topas
wget -nc -O $1/G4Data/G4NDL.4.6.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4NDL.4.6.tar.gz
wget -nc -O $1/G4Data/G4EMLOW.7.9.1.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4EMLOW.7.9.1.tar.gz
wget -nc -O $1/G4Data/G4PhotonEvaporation.5.5.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4PhotonEvaporation.5.5.tar.gz
wget -nc -O $1/G4Data/G4RadioactiveDecay.5.4.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4RadioactiveDecay.5.4.tar.gz
wget -nc -O $1/G4Data/G4PARTICLEXS.2.1.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4PARTICLEXS.2.1.tar.gz
wget -nc -O $1/G4Data/G4SAIDDATA.2.0.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4SAIDDATA.2.0.tar.gz
wget -nc -O $1/G4Data/G4ABLA.3.1.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4ABLA.3.1.tar.gz
wget -nc -O $1/G4Data/G4INCL.1.0.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4INCL.1.0.tar.gz
wget -nc -O $1/G4Data/G4PII.1.3.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4PII.1.3.tar.gz
wget -nc -O $1/G4Data/G4ENSDFSTATE.2.2.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4ENSDFSTATE.2.2.tar.gz
wget -nc -O $1/G4Data/G4RealSurface.2.1.1.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4RealSurface.2.1.1.tar.gz
wget -nc -O $1/G4Data/G4TENDL.1.3.2.tar.gz http://geant4-data.web.cern.ch/geant4-data/datasets/G4TENDL.1.3.2.tar.gz

[ ! -d "$1/G4Data/G4NDL4.6" ]     && tar -zxvf $1/G4Data/G4NDL.4.6.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4EMLOW7.9.1" ] && tar -zxvf $1/G4Data/G4EMLOW.7.9.1.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/PhotonEvaporation5.5" ] && tar -zxvf $1/G4Data/G4PhotonEvaporation.5.5.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/RadioactiveDecay5.4" ] && tar -zxvf $1/G4Data/G4RadioactiveDecay.5.4.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4PARTICLEXS2.1" ] && tar -zxvf $1/G4Data/G4PARTICLEXS.2.1.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4SAIDDATA2.0" ] && tar -zxvf $1/G4Data/G4SAIDDATA.2.0.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4ABLA3.1" ] && tar -zxvf $1/G4Data/G4ABLA.3.1.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4INCL1.0" ] && tar -zxvf $1/G4Data/G4INCL.1.0.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4PII1.3" ] && tar -zxvf $1/G4Data/G4PII.1.3.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4ENSDFSTATE2.2" ] && tar -zxvf $1/G4Data/G4ENSDFSTATE.2.2.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/RealSurface2.1.1" ] && tar -zxvf $1/G4Data/G4RealSurface.2.1.1.tar.gz -C $1/G4Data
[ ! -d "$1/G4Data/G4TENDL1.3.2" ] && tar -zxvf $1/G4Data/G4TENDL.1.3.2.tar.gz -C $1/G4Data
