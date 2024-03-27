function [LETmaskDirty,LETmaskClean,mLET] = matRad_calcLETmask(dij)
% Calculates logical matrix where LET is above (LETmask.Dirty) or below (LETmask.Clean) a certain threshold
%
% call
%   LETmask = matRad_calcLETmask(dij)
%
% input
%   dij:       matRad dij struct
%
% output
%   LETmask:   logical matrix for dirty dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% instead of:
% mLET = full(dij.mLETDose{1}(:,:))./full(dij.physicalDose{1}(:,:));
% LETmask = mLET > dij.dirtyDoseThreshold;
% logical matrix is calculated like this:

    [i,j,v] = find(dij.physicalDose{1});
    idx = sub2ind(size(dij.physicalDose{1}),i,j);

    mLET = full(dij.mLETDose{1}(idx) ./ v);  
    
    subIxDirty = mLET > dij.dirtyDoseThreshold;
    subIxClean = ~subIxDirty;

    mLET = sparse(i,j,mLET,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    LETmaskDirty = sparse(i(subIxDirty),j(subIxDirty),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);
    LETmaskClean = sparse(i(subIxClean),j(subIxClean),true,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);


end