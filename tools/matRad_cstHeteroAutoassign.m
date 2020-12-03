function [cstHetero] = matRad_cstHeteroAutoassign(cst)
% Prepares the cst file for the heterogeneity correction algorithm
%
% call
%   cstHetero = matRad_cstHeteroAutoassign(cst)
%
% input
%   cst:            matRad cst struct
%
% output
%   cstHetero:      updated matRad cst struct with 'Lung' property
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cstHetero = cst;
% automatically assign the 'Lung' property to these segmentations
% lungTissue={'Lung','GTV','PTV','CTV','ITV'};  
lungTissue={'Lung'};  

for i = 1:length(cst(:,1))
    if contains(cst{i,2},lungTissue)
        cstHetero{i,5}.HeterogeneityCorrection = 'Lung';
    end
end

end