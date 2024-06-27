#!/usr/bin/env bash

sudo chmod +x .github/runtests.sh
sudo chmod +x thirdParty/MCsquare/bin/MCsquare_linux

mv matRad/optimization/optimizer/ipopt.m optimization/optimizer/ipopt.m.bak

octave --no-gui --eval "pkg install -forge dicom"
octave --no-gui --eval "pkg install -forge nan"
