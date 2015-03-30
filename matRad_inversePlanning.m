function optResult = matRad_inversePlanning(dij,cst,pln)
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
%   XXX
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

% intial fluence profile = uniform bixel intensities
wInit = ones(dij.totalNumOfBixels,1);

% precalculate hadamard product of sparse matrices
if pln.bioOptimization == true
   dij.doseSkeleton   = spones(dij.dose);
   dij.mSqrtBetaDose2 = 2*dij.mSqrtBetaDose;
end

% define objective function
if pln.bioOptimization == true
    objFunc =  @(x) matRad_bioObjFunc(x,dij,cst);
else 
    objFunc =  @(x) matRad_objFunc(x,dij,cst);
end

% minimize objetive function
optResult = matRad_optimize(objFunc,wInit);

% calc dose and reshape from 1D vector to 2D array
optResult.physicalDose = reshape(dij.dose*optResult.w,dij.dimensions);

if pln.bioOptimization == true
    
    a_x = zeros(size(optResult.physicalDose));
    b_x = zeros(size(optResult.physicalDose));
    sPrescription=0;
    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,8}) = cst{i,9}.alphaX;
            b_x(cst{i,8}) = cst{i,9}.betaX;
            if isequal(cst{i,3},'TARGET') 
                sPrescription=cst{i,4};
            end
        end
    end
    
    optResult.effect = (dij.mAlphaDose*optResult.w+(dij.mSqrtBetaDose*optResult.w).^2);
    optResult.effect = reshape(optResult.effect,dij.dimensions);
    
    optResult.RBEWeightedDose = ((sqrt(a_x.^2 + 4 .* b_x .* optResult.effect) - a_x)./(2.*b_x));
    
    %only consider RBE at certain dose level
    sCutOffPhysDose = 0.05;
    optResult.RBE = optResult.RBEWeightedDose./optResult.physicalDose;
    optResult.RBE(optResult.physicalDose < sPrescription*sCutOffPhysDose)=0;

    % a different way to calculate RBE is as follows - leads to the same
    %optResult.RBE = ((sqrt(a_x.^2 + 4 .* b_x .* optResult.effect) - a_x)./(2.*b_x.*optResult.physicalDose));
    %optResult.RBE= reshape(optResult.RBE,dij.dimensions);
      
    optResult.alpha = (dij.mAlphaDose.*spfun(@(x)1./x,dij.dose)) * optResult.w;
    optResult.alpha = reshape(optResult.alpha,dij.dimensions);
    optResult.beta = ( (dij.mSqrtBetaDose.*spfun(@(x)1./x,dij.dose)) * optResult.w ).^2;
    optResult.beta = reshape(optResult.beta,dij.dimensions);
    
    
end

% Make a sound when finished.
beep;
