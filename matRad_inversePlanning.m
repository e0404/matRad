function [wOpt,dOpt] = matRad_inversePlanning(dij,cst,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call [w,dOpt] = mPlan2D_inversePlanning(dij,voi,cst)
% to optimize the intensity modulation according to the constraints
% specified in cst
% wOpt: optimized bixel weight vector
% dOpt: optimized fluence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate meta information for optimization
%optInfoArrays = mr_genOptInfoArrays(cst);

% intial fluence profile = uniform bixel intensities
wInit = ones(dij.totalNumOfBixels,1);

% define objective function

if pln.bioOptimization == false 
    objFunc =  @(x) matRad_IMRTObjFunc(x,dij.dose,cst);
else
    objFunc =  @(x) matRad_IMRTBioObjFunc(x,dij,cst);
end

%[f_ref, g_ref, ~] = matRad_IMRTObjFunc(wInit,dij.dose,cst);
% [f, g, ~] = matRad_IMRTBioObjFunc(wInit,dij,cst);
% % test gradient
% epsilon = 0.005;
% for i = 1:numel(wInit)
%     
%     wDelta = wInit;
%     wDelta(i) = wDelta(i) + epsilon;
%     [fDelta, ~, ~] = matRad_IMRTBioObjFunc(wDelta,dij,cst);
%     
%     numGrad = (fDelta-f)/epsilon;
%     
%     fprintf(['Component # ' num2str(i) ' - percent diff in numerical and analytical gradient = ' ...
%         num2str((numGrad/g(i)-1)*100) '\n']);   
% end

% minimize objetive function
[wOpt,dOpt] = matRad_optimize(objFunc,wInit);

% reshape from 1D vector to 2D array
dOpt = reshape(dOpt,dij.dimensions);

% Make a sound when finished.
beep;