function optResult = matRad_directApertureOptimization(dij,cst,apertureInfo,optResult,visBool,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to run direct aperture optimization
%
% call
%   optResult = matRad_directApertureOptimization(dij,cst,apertureInfo,visBool,varargin)   
%
% input
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   shaoInfo:   aperture shape info struct
%   optResult:  resultGUI struct to which the output data will be added, if
%               this field is empty optResult struct will be created
%               (optional)
%   visBool:    plots the objective function value in dependence of the
%               number of iterations
%   varargin:   optinal: convergence criteria for optimization and biological
%               optimization mode
%
% output
%   optResult:  struct containing optimized fluence vector, dose, and shape
%               info
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
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

if nargin < 5
    visBool = 0;
end

% adjust overlap priorities
cst = matRad_setOverlapPriorities(cst);

% get initial vector for optimization
initApertureInfoVec = apertureInfo.apertureVector;

% set objective function
objFunc =  @(x) matRad_daoObjFunc(x,apertureInfo,dij,cst);

% set projection function
projFunc = @(x) matRad_daoProjectionFunction(x,apertureInfo.limMx,apertureInfo.totalNumOfShapes);

% verify gradients
%matRad_verifyGradient(objFunc,initapertureInfoVec);

% minimize objetive function
optApertureInfoVec = matRad_projectedLBFGS(objFunc,projFunc,initApertureInfoVec,visBool,varargin);

% update the apertureInfoStruct and calculate bixel weights
optResult.apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,optApertureInfoVec);

% override also bixel weight vector in optResult struct
optResult.w    = optResult.apertureInfo.bixelWeights;
optResult.wDao = optResult.apertureInfo.bixelWeights;
% calc dose and reshape from 1D vector to 3D array
optResult.physicalDose = reshape(dij.physicalDose*optResult.w,dij.dimensions);
