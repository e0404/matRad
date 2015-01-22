# General information

matRad is an open source treatment planning system for radiation therapy
written in Matlab. It is meant for educational purposes and supports 
planning of intensity-modulated radiation therapy for mutliple 
modalities. The source code is maintained by a development team at the 
German Cancer Reserach Center DKFZ in Heidelberg, Germany, and other
contributors around the world. We are always looking for more people
willing to help improve matRad. Do not hesitate and get in touch.

---

Copyright 2015, Mark Bangert, on behalf of the matRad development team

m.bangert@dkfz.de

matrad is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option)
any later version.

matRad is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License in the
file license.txt along with matRad. If not, see
<http://www.gnu.org/licenses/>.

---

We are currently working on the first beta release which is anticiapted
in spring 2015.

---

# Preliminary programming rules
* use a prefix 'matRad_' before all matrad code files!
* use capital letters to indicate a new word within variable names, e.g. numOfBeams not (!) num_of_beams
* use extension '_vox' for every varibale that specifies a location in voxel coordinates!
* all visualization code should be in one location within an %.m file
* if we have multiple visulaizations use a variable visMode = 0,1,2,...
