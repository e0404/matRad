[![Current Release](https://img.shields.io/github/v/release/e0404/matRad)](https://github.com/e0404/matRad/releases) 
[![Downloads](https://img.shields.io/github/downloads/e0404/matRad/total)](https://github.com/e0404/matRad/releases) 
[![Contributors](https://img.shields.io/github/contributors/e0404/matRad)](https://github.com/e0404/matRad/graphs/contributors)
[![GitPitch](https://gitpitch.com/assets/badge.svg)](https://gitpitch.com/e0404/matRad/master) 

[![TravisCI Build Status](https://travis-ci.org/e0404/matRad.svg?branch=dev)](https://travis-ci.org/e0404/matRad)
[![Azure Pipelines Build Status](https://dev.azure.com/e0404/matRad/_apis/build/status/e0404.matRad?branchName=dev)](https://dev.azure.com/e0404/matRad)

DOIs:
- General DOI: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3879615.svg)](https://doi.org/10.5281/zenodo.3879615)
- Latest Release: [![DOI](https://zenodo.org/badge/29671667.svg)](https://zenodo.org/badge/latestdoi/29671667)

# General information

---

matRad is an open source treatment planning system for radiation therapy written in Matlab. It supports planning of intensity-modulated radiation therapy for mutliple modalities and is meant **for educational and research purposes**. **IT IS NOT SUITABLE FOR CLINICAL USE** (also see the no-warranty clause in the GPL license). The source code is maintained by a development team at the German Cancer Reserach Center - DKFZ in Heidelberg, Germany, and other contributors around the world. We are always looking for more people willing to help improve matRad. Do not hesitate and [get in touch](mailto:contact@matRad.org).

More information can be found on the project page  at <http://e0404.github.io/matRad/>; a wiki documentation is under constant development at <https://github.com/e0404/matRad/wiki>.

# Getting Started

It’s the first time you want to use matRad?

First, get a local copy of matRad by download or git cloning. Having done that, we recommend you navigate into the folder in Matlab and execute 
```
matRad_rc
```
which will setup the path & configuration and tell you the current version.

Then there’re three options for a pleasant start with matRad. Choose one or try out each of them.

### Option 1: Using the GUI

For an intuitive workflow with the graphical user interface, type 
```
matRadGUI
```
in your command window. An empty GUI should be opened. Click the _*Load.mat_ data-Button in the Workflow-section to load a patient. Set the plan and optimization parameters, calculate the dose influence matrix and execute the fluence optimization in the GUI.

### Option 2: Using the main script

If you prefer scripting, open the default script *matRad.m* from the main matRad folder:
```
edit matRad.m
```
Use it to learn something about the code structure and execute it section by section.

You can also run the full script for an example photon plan by just typing
```
matRad
``` 
in your command window.

### Option 3: Using the examples

The most time consuming but also most educational approach to matRad. 

When in the main matRad folder, navigate to the folder *examples*. Open one of the examples given there. Execute it section by section. Move on to the next example afterwards.

# Need help?

If you encounter problems with matRad, please consider the following guidelines **before** submitting issues on our github page. 

* Check you are using the newest version of matRad.
* Please check the description of how to set up matRad and its technical documentation in the [wiki](https://github.com/e0404/matRad/wiki).
* Go through the relevant examples and see if they answer your question (see *Option 3* above!)
* Check open and closed issues for your question.

Still having problems? Then create an issue, provide a **minimum example** of your attempted workflow / what causes the problems and be patient!

# Citing matRad

### Scientific papers

If you use matRad in a scientific publication, consider citing the following paper:

Wieser, Hans-Peter, et al. "Development of the open-source dose calculation and optimization toolkit matRad." Medical Physics 44.6 (2017): 2556-2568. 

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fmp.12251-blue)](https://doi.org/10.1002/mp.12251) 

BibTex entry:
```
@article{wieser2017development,
  title={Development of the open-source dose calculation and optimization toolkit matRad},
  author={Wieser, Hans-Peter and Cisternas, Eduardo and Wahl, Niklas and Ulrich, Silke and Stadler, Alexander and Mescher, Henning and M{\"u}ller, Lucas-Raphael and Klinge, Thomas and Gabrys, Hubert and Burigo, Lucas and others},
  journal={Medical Physics},
  volume={44},
  number={6},
  pages={2556--2568},
  year={2017},
  publisher={Wiley Online Library},
  doi={10.1002/mp.12251}
}
```

### Citing as Software

matRad's code also has its own general DOI with Zenodo: 

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3879615.svg)](https://doi.org/10.5281/zenodo.3879615)

You can cite specific versions of matRad in your work! For example, Here is the badge that lead's to the latest release of matRad:

[![DOI](https://zenodo.org/badge/29671667.svg)](https://zenodo.org/badge/latestdoi/29671667)

---

Copyright 2020 the matRad development team. 

matrad@dkfz.de

All the elements of the compilation of matRad and Ipopt are free software. You can redistribute and/or modify matRad's source code version provided as files with .m and .mat extension under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 (GPL v3). You can also add to matRad the Ipopt functionality by using the precompiled mex files of the Ipopt optimizer in object code version which are licensed under the Eclipse Public License Version 1.0 (EPL v1.0), also made available for download via https://projects.coin-or.org/Ipopt.
matRad also contains interfaces to an open-source photon Monte Carlo dose calculation engine developed by Edgardo Dörner hosted on GitHub (http://github.com/edoerner/ompMC) and to the open-source proton Monte Carlo project MCsquare (www.openmcsquare.org) from UCLouvain, Louvain-la-Neuve, Belgium. Both interfaces are integrated into matRad as submodules.

In addition, we provide a matlab standalone version of the compilation of matRad and Ipopt, where the files of matRad and Ipopt are licensed under GPL v3 and EPL v1.0 respectively. The matlab standalone version is meant to be used by students for learning and practicing scientific programming and does not yet contain the interfaces to the aforementioned Monte Carlo dose calculation engines.

matRad is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Please note that we treat the compilation of matRad and Ipopt as separate and independent works (or modules, components, programs). Therefore, to the best of our understanding, the compilation of matRad and Ipopt is subject to the "Mere Aggregation" exception in section 5 of the GNU v3 and the exemption from "Contributions" in section 1. b) ii) of the EPL v1.0. Should this interpretation turn out to be not in compliance with the applicable laws in force, we have provided you with an additional permission under GNU GPL version 3 section 7 to allow you to use the work resulting from combining matRad with Ipopt.

You will receive a copy of the GPL v3  and a copy of the EPL v1.0 in the file LICENSES.txt along with the compilation. If not, see http://www.gnu.org/licenses/ and/or http://opensource.org/licenses/EPL-1.0/.

---
