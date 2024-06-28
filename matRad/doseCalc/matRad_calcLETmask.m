function [LETmaskDirty,LETmaskClean,mLET] = matRad_calcLETmask(dij)
% Calculates logical matrix where LET is above (LETmaskDirty) or below
% (LETmaskClean) a certain threshold which is set in matRad_calcDirtyDose
%
% call
%   [LETmaskDirty,LETmaskClean,mLET] = matRad_calcLETmask(dij)
%
% input
%   dij:       matRad dij struct after the dose calculation using
%              matRad_calcParticleDose
%
% output
%   LETmaskDirty:   logical matrix (0 and 1) for dirty dose 
%                   -> it needs to go into matRad_calcDirtyDose to create
%                   dirty dose as a variable
%   LETmaskClean:   logical matrix (0 and 1) for clean dose -> into matRad_calcDirtyDose
%   mLET:           mLETDose matrix divided by the physicalDose matrix
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try %Leverage similar Nonzero Structure    
    for scenIx = find(~cellfun(@isempty,dij.physicalDose))
        [i,j] = find(dij.mLETDose{scenIx});
        mLET{scenIx} = sparse(i,j,nonzeros(dij.mLETDose{scenIx}) ./ nonzeros(dij.physicalDose{scenIx}),dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    
        LETmaskDirty{scenIx} = mLET{scenIx} > dij.dirtyDoseThreshold;
        subIxClean = nonzeros(mLET{scenIx}) < dij.dirtyDoseThreshold; 
        
        LETmaskClean{scenIx} = sparse(i(subIxClean),j(subIxClean),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
        
        %{
        [i,j,v] = find(dij.physicalDose{1});
        idx = sub2ind(size(dij.physicalDose{1}),i,j);

        mLET{1} = full(dij.mLETDose{1}(idx) ./ v);

        subIxDirty = mLET{1} > dij.dirtyDoseThreshold;
        subIxClean = ~subIxDirty;

        mLET{1} = sparse(i,j, mLET{1},dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);

        LETmaskDirty{1} = sparse(i(subIxDirty),j(subIxDirty),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
        LETmaskClean{1} = sparse(i(subIxClean),j(subIxClean),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
        %}
    end
catch ME
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Nonzero structure of LETxDose and physicalDose matrix not similar or matrices do not exist!');    
end

    
end