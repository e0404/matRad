% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to compile mex files. All mex files will be stored in
% matRadRootDir/APM/MEX
% 
% call
%   matRad_compileAPMmexFiles();
%
% input
%  
% output
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 Hans-Peter Wieser
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clc,clear,close 

cfg = coder.config('mex');
% if isunix
%     cfg.PostCodeGenCommand = 'buildInfo.addLinkFlags(''-fopenmp'')'; 
% end

codegen -config cfg matRad_sumGauss 

codegen -config cfg matRad_calcGeoDistsFast -report 

codegen -config cfg matRad_calcSecLatMom -report 
codegen -config cfg matRad_calcSecRangeMom -report

codegen -config cfg matRad_calcSecLatMomFast -report 
codegen -config cfg matRad_calcSecRangeMomFast -report

codegen -config cfg matRad_calcThirdRangeMom -report 
codegen -config cfg matRad_calcFourthRangeMom -report





