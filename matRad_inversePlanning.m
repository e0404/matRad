function [wOpt,dOpt] = matRad_inversePlanning(dij,cst,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
% 
% call
%   [wOpt,dOpt] = matRad_inversePlanning(dij,cst)
%
% input
%   dij:    matRad dij struct
%   cst:    matRad cst struct
%   pln:    matRad pln struct
%
% output
%   wOpt:   optimized bixel weight vector
%   dOpt:   optimized dose distribution (as cube)
%
% References
%   -
%
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

%precalculate hadamard product of sparse matrices
if pln.bioOptimization == true && strcmp(pln.radiationMode,'carbon')
   dij.doseSkeleton = spones(dij.dose);
   dij.mAlphaDose = dij.mAlpha.*dij.dose;
   dij.mBetaDose = sqrt(dij.mBeta).*dij.dose;
end
% define objective function
if pln.bioOptimization == true && strcmp(pln.radiationMode,'carbon')
    objFunc =  @(x) matRad_IMRTBioObjFunc(x,dij,cst);
else 
    objFunc =  @(x) matRad_IMRTObjFunc(x,dij.dose,cst);
end

% minimize objetive function
[wOpt,dOpt] = matRad_optimize(objFunc,wInit);

% reshape from 1D vector to 2D array
dOpt.PhysicalDose = reshape(dOpt.PhysicalDose,dij.dimensions);

if pln.bioOptimization == true && strcmp(pln.radiationMode,'carbon')    
    a_x = zeros(size(dOpt.Effect,1),1);
    b_x = zeros(size(dOpt.Effect,1),1);
    
    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            ind =find([dij.baseData(1).Tissue.Class]==cst{i,9}.TissueClass);
            a_x(cst{i,8})=dij.baseData(1).Tissue(ind).alphaX;
            b_x(cst{i,8})=dij.baseData(1).Tissue(ind).betaX;
        end
    end
    
    dOpt.BiologicalDose = ((sqrt(a_x.^2 + 4 .* b_x .* dOpt.Effect) - a_x)./(2.*b_x));
    dOpt.BiologicalDose = reshape(dOpt.BiologicalDose,dij.dimensions);
    
    dOpt.RBE = dOpt.BiologicalDose./dOpt.PhysicalDose;
    % a different way to calculate RBE is as follows - leads to the same
    %dOpt.RBE = ((sqrt(a_x.^2 + 4 .* b_x .* dOpt.Effect') - a_x)./(2.*b_x.*dOpt.PhysicalDose'))';
    %dOpt.RBE= reshape(dOpt.RBE,dij.dimensions);
    dOpt.Effect = reshape(dOpt.Effect,dij.dimensions);
    dOpt.Alpha = reshape(dij.mAlpha*wOpt,dij.dimensions);
end


% Make a sound when finished.
beep;