function [dij] = matRad_calculateProbabilisticQuantities(dij,cst,pln,mode4D)
% matRad helper function to compute probabilistic quantities, i.e. expected
%   dose influence and variance omega matrices for probabilistic
%   optimization
%
% call:
%   [dij,cst] = matRad_calculateProbabilisticQuantities(dij,cst,pln)
%
% input:
%   dij:        matRad dij struct
%   cst:        matRad cst struct (in dose grid resolution)
%   pln:        matRad pln struct
%   mode4D:     (Optional) Handling of 4D phases: 
%               - 'all'   : include 4D scen in statistic
%               - 'phase' : create statistics per phase (default)
%
% output:
%   dij:        dij with added probabilistic quantities
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
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

if ~pln.bioModel.bioOpt
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

% Full scenario indices in the order of the scenario numbers, so that they
% match pln.multScen.scenProb and pln.multScen.linearMask
scens = pln.multScen.scenNum2scenIx();
scenProb = pln.multScen.scenProb;

%Create structures
for i = 1:numel(fNames)
    matRad_cfg.dispInfo('\tE[D] & Omega[D] for ''%s'':\n',fNames{1,i});
    if strcmpi(mode4D,'phase')
        % expected quantities are indexed with the CT scenario id
        numOfCtIx = size(pln.multScen.scenMask,1);
        ctIxMap = pln.multScen.linearMask(:,1);
    else
        numOfCtIx = 1;
        ctIxMap = ones(numel(scens),1);
    end
    usedCtIx = unique(ctIxMap)';

    dij.([fNames{1,i} 'Exp']) = cell(1,numOfCtIx);
    dij.([fNames{1,i} 'Exp'])(usedCtIx) = {spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1)};
    dij.([fNames{1,i} 'Omega']) = cell(size(cst,1),numOfCtIx);
    [dij.([fNames{1,i} 'Omega']){voiIx,usedCtIx}] = deal(zeros(dij.totalNumOfBixels));

    %Now loop over the scenarios
    
    for s = 1:numel(scens)
        matRad_cfg.dispInfo('\t\tScenario %d/%d...',s,numel(scens));
        ctIx = ctIxMap(s);
        scenIx = scens(s);
        
        %Add up Expected value
        dij.([fNames{1,i} 'Exp']){ctIx} = dij.([fNames{1,i} 'Exp']){ctIx} + dij.([fNames{1,i}]){scenIx} .* scenProb(s);

        %Add up Omega
         for v = voiIx
             dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} + ...
                 (((dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)' * scenProb(s)) * ...
                 (dij.(fNames{1,i}){scenIx}(cst{v,4}{ctIx},:)) * scenProb(s)));
         end
         matRad_cfg.dispInfo('Done!\n');
    end
    matRad_cfg.dispInfo('\tFinalizing Omega...');
    %Finalize Omega matrices
    for ctIx = usedCtIx
        for v = voiIx
            dij.([fNames{1,i} 'Omega']){v,ctIx} = dij.([fNames{1,i} 'Omega']){v,ctIx} - (dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:)' * dij.([fNames{1,i} 'Exp']){ctIx}(cst{v,4}{ctIx},:));
        end
    end  
    matRad_cfg.dispInfo('\tDone!\n');
end
