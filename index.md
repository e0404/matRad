---
layout: default
---

![matRad](/assets/img/matrad_hat.svg){:height="10%" width="10%"} is an open source software for radiation treatment planning of intensity-modulated photon, proton, and carbon ion therapy. matRad is developed for educational and research purposes. It is entirely written in [MATLAB](http://www.mathworks.com/products/matlab).

For more information about this project check out the [matRad overview](https://github.com/e0404/matRad/wiki/documents/matRad.pdf) or watch a recodring of the matRad webinar at the brown bag medical physics seminar at Massachusetts General Hospital.

<iframe width="100%" height="300px" src="https://www.youtube.com/embed/40_n7BIqLdw" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
&nbsp;

Help us fund matRad
========

The sustainable development of matRad costs time. To help us make a strong case towards our funders that we should continue to allocate resources for matRad's development in the future, please <a id="writeEmail" href="mailto:contact@matrad.org?Subject=We use matRad at ..." target="_top">do email</a> us and let us know if you are using matRad in your department. We maintain a list of current user groups in a <a id="gotoMap" href="https://drive.google.com/open?id=1nG8_pE0RR5j1lp2iLLeBCeuwrO0&usp=sharing">custom google map</a>. This helps us to show the relevance of our work for the radiation oncology community. As we do not get detailed download statistics from Github, such information is very valuable to us. Thank you!

Also, if you are interested in irregular updates, please send us an <a href="mailto:contact@matrad.org?Subject=matRad newsletter" target="_top">email</a> to receive the matRad newsletter.

matRad developments were/are, in part, supported by the following research grants:
- Grant  No. BA 2279/3-1 (Project No. 265744405) from the German Research Foundation (DFG)
- Grant No. WA 4707/1-1 (Project No. 443188743) from the German Research Foundation (DFG)
- Grant No. 70113094 from the German Cancer Aid


Features
========

matRad comprises

* MATLAB functions to model the entire treatment planning workflow
* Example patient data
* Physical and biological base data for all required computations

In particular we provide functionalities for

* Ray tracing
* Photon dose calculation
* Proton dose calculation
* Carbon ion dose calculation (including 3D RBE modeling)
* Inverse planning (based on physical dose and biological effect)
* Multileaf collimator sequencing
* Basic treatment plan visualization and evaluation

matRad is constantly evolving. If you are looking for a special feature do not hesitate and <a href="mailto:contact@matrad.org?Subject=matRad" target="_top">get in touch</a>.

Documentation
========

A [Wiki](https://github.com/e0404/matRad/wiki) documentation is under constant development.

Contact
========
If you have any questions or wish to contribute to the development of matRad you can either write an email to <a href="mailto:contact@matrad.org">contact@matrad.org</a>, use the GitHub functionalities to file a new issue, or directly fork the matRad repository and send pull requests with your own custom developments. What are you waiting for?

<a class="twitter-timeline" href="https://twitter.com/mat_rad?ref_src=twsrc%5Etfw">Tweets by mat_rad</a> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>


Development team
========
matRad development is driven by the research group [Radiotherapy Optimization](https://www.dkfz.de/en/medphys/optimization_algorithms/optimization_algorithms.html) within the [Division of Medical Physics in Radiation Oncology](https://www.dkfz.de/en/medphys/) at the [German Cancer Research Center DKFZ](https://www.dkfz.de/en/index.html) in [Heidelberg](https://www.google.de/maps/place/German+Cancer+Research+Center/@49.41023,8.6730541,5129m/data=!3m1!1e3!4m5!3m4!1s0x4797c12f081fefc5:0x21c080e0d5b438c8!8m2!3d49.4144184!4d8.6719509?hl=en). Please see our [author list](https://github.com/e0404/matRad/blob/master/AUTHORS.txt) for an overview of all code contributers; also beyond there are many people supporting matRad through consulting and provision of data.

We are actively looking for beta testers that can provide external feedback on our code and developers that would like to take an active role in the future. Do not hesitate and <a href="mailto:contact@matrad.org?Subject=matRad" target="_top">get in touch</a>.

Scientific publications
========
matRad has been used for the conduct of several [peer-reviewed journal publications](https://scholar.google.de/scholar?cites=6882863808740897521).

License and disclaimer
========
matRad is offered as a compilation with the Ipopt software package for large-scale nonlinear optimization. All the elements of the compilation of matRad and Ipopt are free software. You can redistribute and/or modify matRad's source code provided as files with .m and .mat extension under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 (GPL v3). You can also add to matRad the Ipopt functionality by using the precompiled mex files of the Ipopt optimizer in object code version which are licensed under the Eclipse Public License Version 1.0 (EPL v1.0), also made available for download via [https://github.com/coin-or/Ipopt](https://github.com/coin-or/Ipopt).

In addition, we provide a matlab standalone version of the compilation of matRad and Ipopt, where the files of matRad and Ipopt are licensed under GPL v3 and EPL v1.0 respectively. The matlab standalone version is meant to be used by students for learning and practicing scientific programming.

matRad is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Please note that we treat the compilation of matRad and Ipopt as separate and independent works (or modules, components, programs). Therefore, to the best of our understanding, the compilation of matRad and Ipopt is subject to the "Mere Aggregation" exception in section 5 of the GNU v3 and the exemption from "Contributions" in section 1. b) ii) of the EPL v1.0. Should this interpretation turn out to be not in compliance with the applicable laws in force, we have provided you with an additional permission under GNU GPL version 3 section 7 to allow you to use the work resulting from combining matRad with Ipopt.

You will receive a copy of the GPL v3 and a copy of the EPL v1.0 in the file LICENSES.txt along with the compilation. If not, see [http://www.gnu.org/licenses/](https://www.gnu.org/licenses/licenses.en.html) and/or [http://opensource.org/licenses/EPL-1.0/](http://opensource.org/licenses/EPL-1.0/).







