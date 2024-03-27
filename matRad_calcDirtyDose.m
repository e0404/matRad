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

if ~exist('dij','var') || isempty(dij)
    disp('Not enough input arguments! Calculation is not working.')
else
    if ~exist('LET_thres','var') || isempty(LET_thres)
        disp('Not enough input arguments! Calculation is not working.')
    else
        dij.dirtyDoseThreshold                       = LET_thres;
        [dij.LETmaskDirty,dij.LETmaskClean,dij.mLET] = matRad_calcLETmask(dij);
        dij.dirtyDose{1}                             = dij.LETmaskDirty .* dij.physicalDose{1};
        dij.cleanDose{1}                             = dij.LETmaskClean .* dij.physicalDose{1};
    end
end

end