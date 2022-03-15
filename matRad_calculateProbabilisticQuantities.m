function [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln,mode4D)
% matRad helper function to compute probabilistic quantities, i.e. expected
%   dose influence and variance omega matrices for probabilistic
%   optimization
%
% call
%   [dij,cst] = matRad_calculateProbabilisticQuantities(dij,cst,pln)
%
% input
%   dij:        matRad dij struct
%   cst:        matRad cst struct (in dose grid resolution)
%   pln:        matRad pln struct
%   mode4D:     (Optional) Handling of 4D phases: 
%               - 'all'   : include 4D scen in statistic
%               - 'phase' : create statistics per phase (default)
%
% output
%   dij:        dij with added probabilistic quantities
%
% References
%   -
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

matRad_cfg = MatRad_Config.instance();

if nargin < 4
    mode4D = 'phase';
elseif ~any(strcmpi(mode4D,{'phase','all'}))
    matRad_cfg.dispError('4D calculation mode %s unknown, allowed is only ''phase'' or ''all''',mode4D);
end
    

matRad_cfg.dispInfo('Calculating probabilistic quantities E[D] & Omega[D] ...\n');

if ~pln.bioParam.bioOpt
    fNames = {'physicalDose'};
else
    fNames = {'mAlphaDose','mSqrtBetaDose'};
end

% create placeholders for expected influence matrices

voiIx = [];
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        voiIx = [voiIx i];
    end
end

scens = find(pln.multScen.scenMask);

%Create structures
for i = 1:numel(fNames)
    matRad_cfg.dispInfo('\tE[D] & Omega[D] for ''%s'':\n',fNames{1,i});
    if strcmpi(mode4D,'phase')
        dij.([fNames{1,i} 'Exp']){1:pln.multScen.numOfCtScen} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
        dij.([fNames{1,i} 'Omega']) = cell(size(cst,1),pln.multScen.numOfCtScen);
        
        ixTmp = cell(ndims(pln.multScen.scenMask),1);
        [ixTmp{:}] = ind2sub(size(pln.multScen.scenMask),scens);
        ctIxMap = ixTmp{1};        
    else
        dij.([fNames{1,i} 'Exp']){1} = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);
        dij.([fNames{1,i} 'Omega']) = cell(size(cst,1),1);
        ctIxMap = ones(numel(scens),1);
    end
    [dij.([fNames{1,i} 'Omega']){voiIx,:}] = deal(zeros(dij.totalNumOfBixels));
    
    %Now loop over the scenarios
    
    for s = 1:numel(scens)
        matRad_cfg.dispInfo('\t\tScenario %d/%d...',s,numel(scens));
        ctIx = ctIxMap(s);
        scenIx = scens(s);
        
        %Add up Expected value
        dij.([fNames{1,i} 'Exp']){ctIx} = dij.([fNames{ctIx,i} 'Exp']){ctIx} + dij.([fNames{1,i}]){scenIx} .* pln.multScen.scenProb(s);
        
        %Add up Omega
         for v = voiIx
             dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} + ...
                 (((dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)' * pln.multScen.scenProb(s)) * ...
                 (dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)) * pln.multScen.scenProb(s)));
         end
         matRad_cfg.dispInfo('done!\n');
    end
    matRad_cfg.dispInfo('\tFinalizing Omega...');
    %Finalize Omega matrices
    unCtIx = unique(ctIx);    
    for ctIx = unCtIx
        for v = voiIx
            dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} - (dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:)' * dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:));
        end
    end  
    matRad_cfg.dispInfo('\tDone!\n');
end


