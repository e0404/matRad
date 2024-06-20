#!/usr/bin/env bash

sudo chmod +x matRad/.github/runtests.sh
sudo chmod +x matRad/MCsquare/bin/MCsquare_linux

mv optimization/optimizer/ipopt.m optimization/optimizer/ipopt.m.bak

octave --no-gui --eval "pkg install -forge dicom"
octave --no-gui --eval "pkg install -forge nan"
