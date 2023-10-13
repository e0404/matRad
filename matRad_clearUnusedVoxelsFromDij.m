function [dij, mask] = matRad_clearUnusedVoxelsFromDij(cstOnDoseGrid, dij, scenarios)
% matRad function to set the voxels in dij that are not used for
% optimization.
%
% call
%   [dij] = matRad_clearUnusedVoxelsFromDij(cst, dij)
%   [dij] = matRad_clearUnusedVoxelsFromDij(cst, dij, scenarios)
%   [dij, mask] = matRad_clearUnusedVoxelsFromDij(cst, dij)
%   [dij, mask] = matRad_clearUnusedVoxelsFromDij(cst, dij, scenarios)
%
% input
%   cstInDoseGrid:          cst (on dose grid)
%   dij:                    dij struct
%   scenarios (optional):   explicitly define the scenario indexes that need to be cleared
%                               
% 
% output
% 
%   dij:                cleared dij struct
%   mask (optional):    cell array containig the mask that has been applied to every ct scenario
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
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
    matRad_cfg.dispInfo('Clearing unused voxels in dij... ');
    
    % If no scenarios are specified, clear them all
    if ~exist('scenarios', 'var')
        scenarios = find(~cellfun(@isempty, dij.physicalDose));
    end

    %convert selected scenario indexes to subscripts
    scenIndexes = cell(ndims(dij.physicalDose),1);
    [scenIndexes{:}] = ind2sub(size(dij.physicalDose),scenarios);


    % get voxel mask with voxels on those structures that have objectives/constraints
    includeMask = matRad_selectVoxelsFromCst(cstOnDoseGrid,dij.doseGrid,'objectivesOnly');

    ctScenarios = unique(scenIndexes{1});

    for ctIdx=ctScenarios'
    
        idx = find(scenIndexes{1}==ctIdx);

        %get linear indexes of the scenarios in the current ct phase
        scenInPhase = sub2ind(size(dij.physicalDose), scenIndexes{1}(idx), scenIndexes{2}(idx), scenIndexes{3}(idx));
        
        %clear the dij voxels
        for scenIdx=scenInPhase'

            dij.physicalDose{scenIdx} = dij.physicalDose{scenIdx}.*includeMask{ctIdx}(:);
        end

        
        if isfield(dij, 'mAlphaDose')
            dij.mAlphaDose{scenIdx} = dij.mAlphaDose{scenIdx}.*includeMask{ctIdx};
        end
        
        if isfield(dij, 'mSqrtBetaDose')
            dij.mSqrtBetaDose{scenIdx} = dij.mSqrtBetaDose{scenIdx}.*includeMask{ctIdx};
        end

        if isfield(dij,'mLETDose')
            dij.mSqrtBetaDose{scenIdx} = dij.mLETDose{scenIdx}.*includeMask{ctIdx};
        end
    end

    if nargout>1
        mask = includeMask;
    end
    matRad_cfg.dispInfo('done.\n');
end