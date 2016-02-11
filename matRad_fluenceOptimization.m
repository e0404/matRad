function resultGUI = matRad_fluenceOptimization(dij,cst,pln,visBool,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad inverse planning wrapper function
% 
% call
%   resultGUI = matRad_fluenceOptimization(dij,cst,pln,varargin)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct
%   pln:        matRad pln struct
%   visBool:    plots the objective function value in dependence of the
%               number of iterations
%   varargin:   optinal: convergence criteria for optimization and biological
%               optimization mode
%
% output
%   resultGUI struct containing optimized fluence vector, dose, and (for
%   biological optimization) RBE-weighted dose etc.
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

if nargin < 4
    visBool = 0;
end

% intial fluence profile = uniform bixel intensities
wInit = ones(dij.totalNumOfBixels,1);

% consider VOI priorities
cst  = matRad_setOverlapPriorities(cst);

%% adjust internally for fractionation 
for i = 1:size(cst,1)
    for j = 1:size(cst{i,6},2)
        if ~strcmp(cst{i,6}(j).type,'mean') && ~strcmp(cst{i,6}(j).type,'EUD')
          cst{i,6}(j).parameter(2) = cst{i,6}(j).parameter(2)/pln.numOfFractions;
        end
    end
end

% define objective function
if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
        && strcmp(pln.radiationMode,'carbon')
    % check if you are running a supported rad
    dij.ax  = zeros(dij.numOfVoxels,1);
    dij.bx  = zeros(dij.numOfVoxels,1);
    dij.Smax=zeros(dij.numOfVoxels,1);
    
    for i = 1:size(cst,1)
        for j = 1:size(cst{i,6},2)
            % check if only allowed objectives were defined
            if sum(strcmp(cst{i,6}(j).type,{'square overdosing', ...
                                            'square underdosing', ...
                                            'square deviation',...
                                            'mean',...
                                            'EUD'})) < 1
                                        
                error([cst{i,6}(j).type ' objective not supported ' ...
                    'during biological optimization for carbon ions']);
            end
            if strcmp(cst{i,6}(j).type,'square underdosing') || ...
               strcmp(cst{i,6}(j).type,'square deviation')
                if cst{i,6}(j).parameter(2) > 30
                    warning('Prescribed fraction dose > 30Gy. Biological optimization outside the valid domain of the base data.');
                end
            end
        end
        
         if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
             dij.ax(cst{i,4}) = cst{i,5}.alphaX;
             dij.bx(cst{i,4}) = cst{i,5}.betaX;
         end
    end
    
    % define if RBE x dose or biological effect should be used   
    if length(varargin)>1
        IdentifierBioOpt = varargin{1,2};
    else
        IdentifierBioOpt = pln.bioOptimization;
    end
    
    switch IdentifierBioOpt
        case 'effect'
            objFunc = @(x) matRad_bioObjFunc(x,dij,cst);
        case 'RBExD'
            dij.gamma = zeros(dij.numOfVoxels,1);
            idx = dij.bx~=0;  % find indices
            dij.gamma(idx)=dij.ax(idx)./(2*dij.bx(idx)); 
            objFunc = @(x) matRad_bioObjFuncRBExD(x,dij,cst);
    end
  
else
    % set objective function
    objFunc =  @(x) matRad_objFunc(x,dij,cst);
    
end


%% verify gradients
%matRad_verifyGradient(objFunc,dij.totalNumOfBixels);

% set function for projection to feasible set during LBFGS optimization
projFunc = @(x) deal(x.* double(x>0) ,x<=0);

% minimize objective function
wOpt = matRad_projectedLBFGS(objFunc,projFunc,wInit,visBool,varargin);

resultGUI.w = wOpt;
resultGUI.wUnsequenced = wOpt;

% calc dose and reshape from 1D vector to 2D array
resultGUI.physicalDose = reshape(dij.physicalDose*resultGUI.w,dij.dimensions);

if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                            && strcmp(pln.radiationMode,'carbon')

    fprintf('Calculating alpha/beta/effect/cube...');

    a_x = zeros(size(resultGUI.physicalDose));
    b_x = zeros(size(resultGUI.physicalDose));
    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}) = cst{i,5}.alphaX;
            b_x(cst{i,4}) = cst{i,5}.betaX;
        end
    end
    
    resultGUI.effect = (dij.mAlphaDose*resultGUI.w+(dij.mSqrtBetaDose*resultGUI.w).^2);
    resultGUI.effect = reshape(resultGUI.effect,dij.dimensions);
    
    resultGUI.RBExDose = zeros(size(resultGUI.effect));
    ix = resultGUI.effect>0;
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));
    resultGUI.RBE = resultGUI.RBExDose./resultGUI.physicalDose;
   
    AlphaDoseCube    = dij.mAlphaDose * resultGUI.w;
    resultGUI.alpha  = (reshape(AlphaDoseCube,dij.dimensions))./resultGUI.physicalDose;
    SqrtBetaDoseCube = dij.mSqrtBetaDose * resultGUI.w;
    resultGUI.beta   = ((reshape(SqrtBetaDoseCube,dij.dimensions))./resultGUI.physicalDose).^2;
    
    fprintf(' done!\n');
    
end