function dij = matRad_clearUnusedVoxelsFromDij(cst, dij, scenarios)
    
    % If no scenarios are specified, clear them all
    if ~exist('scenarios', 'var')
        scenarios = find(~cellfun(@isempty, dij.physicalDose));
    end

    %convert selected scenario indexes to subscripts
    scenIndexes = cell(ndims(dij.physicalDose),1);
    [scenIndexes{:}] = ind2sub(size(dij.physicalDose),scenarios);

    % get voxel mask with voxels on those structures that have
    % objectives/constraints
    includeMask = matRad_getRobustVoxelsOnGrid(cst,dij.doseGrid,'objectivesOnly');


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
    end
end