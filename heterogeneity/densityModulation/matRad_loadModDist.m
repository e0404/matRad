function [DensMod] = matRad_loadModDist(Pmod)
% matRad function to load density modulation tables
% (only for Pmod = 250, 450 and 800)
%
% call
%   [DensMod] = matRad_loadModDist(Pmod)
%
% input
%   Pmod:           Modulation Power
% 
% output
%   DensMod:        Density probabilty distribution
%
% References
%   [1] https://iopscience.iop.org/article/10.1088/1361-6560/aa641f
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

DensMod = load(['DensMod_',num2str(Pmod),'mu_0,2603_rho_1,5mm_pixelsize.txt']);
DensMod(:,2) = DensMod(:,2) / sum(DensMod(:,2));

end

