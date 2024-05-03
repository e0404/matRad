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

    [i,j,v] = find(dij.physicalDose{1});
    idx = sub2ind(size(dij.physicalDose{1}),i,j);

    mLET = full(dij.mLETDose{1}(idx) ./ v);  
    
    subIxDirty = mLET > dij.dirtyDoseThreshold;
    subIxClean = ~subIxDirty;

    mLET = sparse(i,j,mLET,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    LETmaskDirty = sparse(i(subIxDirty),j(subIxDirty),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    LETmaskClean = sparse(i(subIxClean),j(subIxClean),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);


end