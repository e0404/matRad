#!/usr/bin/env bash

sudo chmod +x .github/runtests.sh
sudo chmod +x thirdParty/MCsquare/bin/MCsquare_linux

mv matRad/optimization/optimizer/ipopt.m optimization/optimizer/ipopt.m.bak

octave --no-gui --eval "pkg install -forge dicom"
octave --no-gui --eval "pkg install -forge nan"
octave --no-gui --eval "pkg install \"https://downloads.sourceforge.net/project/octave/Octave%20Forge%20Packages/Individual%20Package%20Releases/image-2.14.0.tar.gz\""
