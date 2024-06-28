function dij = matRad_calcDirtyDose(LET_thres,dij)
% Calculates Dirty and Clean Dose by using LET threshold
%
% call
%   dij = matRad_calcDirtyDose(LET_thres,dij)
%
% input
%   LET_thres:  LET threshold, above: dirty dose, below: clean dose
%   dij:        matRad dij struct
%
% output
%   dij:        matRad dij struct with dirty and clean dose
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dij.dirtyDoseThreshold                       = LET_thres;
[dij.LETmaskDirty,dij.LETmaskClean,dij.mLET] = matRad_calcLETmask(dij);

dij.dirtyDose = cellfun(@times,dij.LETmaskDirty, dij.physicalDose,'UniformOutput',false);
dij.cleanDose = cellfun(@times,dij.LETmaskClean, dij.physicalDose,'UniformOutput',false);

end