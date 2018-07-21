function [cst, resultGUI] = matRad_calcVar(ct, cst, stf, pln, dij, resultGUI, param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad variance calculation function
%
% call
%   [cst, resultGUI] = matRad_calcVar(ct, cst, stf, pln, dij, resultGUI)
%
% input
%   resultGUI:        matRads resultGUI struct
%   pln:              matRad plan meta information struct
%   dij:              matRad dij dose influence struct
%   stf:              matRad steering information struct
%   cst:              matRad critical structure struct cst
%
% output
%   cst:              matRad critical structure struct cst
%   resultGUI:        matRads resultGUI struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('param','var')
    if ~isfield(param,'logLevel')
        param.logLevel = 1;
    end
else
    param.logLevel    = 1;
end


%% set covariance matrices
pln.multScen.setCovarianceMatrix(ct,cst,pln,stf);

options     = matRad_probOptions(ct,cst,pln);
pln.probOpt = options.probOpt;

if isfield(resultGUI,'wRob')
    param.CALC_OMEGA = false;
    fieldName        = 'physicalDoseExpRob';
    if strcmp(pln.probOpt.InputUCT,'phys')
        % disable covariance matrices to model biological uncertainties
        pln.multScen.mCovBio = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        
    elseif strcmp(pln.probOpt.InputUCT,'bio')
        % disable covariance matrices to model physical uncertainties
        pln.multScen.mCovRangeSys = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovRangeRnd = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovLatSys   = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovLatRnd   = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
    end
    
    
else
    param.CALC_OMEGA = true;
    fieldName        = 'physicalDoseExp';
    
    if strcmp(pln.probOpt.InputUCT,'phys')
        % disable covariance matrices to model biological uncertainties
        pln.multScen.mCovBio = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        
    elseif strcmp(pln.probOpt.InputUCT,'bio')
        % disable covariance matrices to model physical uncertainties
        pln.multScen.mCovRangeSys = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovRangeRnd = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovLatSys   = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
        pln.multScen.mCovLatRnd   = spalloc(dij.totalNumOfBixels,dij.totalNumOfBixels,1);
    end
end

dose = resultGUI.(fieldName)(pln.probOpt.voxelList);

% determine voxel indices for variance calculation
if pln.probOpt.relDoseThreshold > 0
    ix = pln.probOpt.voxelList(dose > pln.probOpt.relDoseThreshold * max(resultGUI.(fieldName)(:)))';
else
    ix = pln.probOpt.voxelList';
end


switch pln.multScen.TYPE
    
    case 'apmScen'
        
        if strcmp(pln.probOpt.typeVarCalc,'serial')
            
            [cst,resultGUI] = matRad_calcParticleVarSerial(ct,cst,stf,pln,dij,resultGUI,ix,param);
            
        elseif strcmp(pln.probOpt.typeVarCalc,'parallel')
            
            [cst,resultGUI] = matRad_calcParticleVarParallel(ct,cst,stf,pln,dij,resultGUI,ix,param);
            
        elseif strcmp(pln.probOpt.typeVarCalc,'parallelParFor')
            
            [cst,resultGUI] = matRad_calcParticleVarParallelParFor(ct,cst,stf,pln,dij,resultGUI,ix,param);
            
        else
            
            matRad_dispToConsole('matRad_calcVar: This option is not implemented!',[],'error')
            
        end
        
    case  {'wcScen','impScen'}
        
        [cst,resultGUI] = matRad_estimateParticleVar(cst,pln,dij,resultGUI);
        
    otherwise
        
        % create placeholde
end






end









