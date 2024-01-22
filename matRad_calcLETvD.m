function dij = matRad_calcLETvD(dij)
% calculates matrix by dividing mLETDose matrix by the physicalDose matrix
[i,j,v] = find(dij.physicalDose{1});
idx = sub2ind(size(dij.physicalDose{1}),i,j);

LETvD = full(dij.mLETDose{1}(idx) ./ v);
LETvD = sparse(i,j,LETvD,dij.doseGrid.numOfVoxels,dij.totalNumOfBixels);

dij.LETvD = LETvD;

end