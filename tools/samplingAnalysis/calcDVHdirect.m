function [dvhPoints, dvhV] =  calcDVHdirect(cst, doseCube, DVHType, dvhPoints)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ddirect dvh calculation
% 
% call
%   calcDVHdirect(cst, doseCube, DVHType, dvhPoints)
%
% input
%   cst:             matRad cst struct
%   doseCube:        matRad dose cube
%   DVHType:         type 1: cumulative DVH
%                    type 2: differential
%   dvhPoints        dose values for DVH. Provide dvhPoints from nominal
%                    scenario to compare multiple dvh
%
% output
%   numerical calculation of DVH without graphical output   
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DVH cummulative (1) / differential (2)
if ~exist('DVHType','var') || isempty(DVHType)
    DVHType = 1;
end

if ~exist('dvhPoints', 'var')
    maxDose = max(doseCube(:));
    minDose = min(doseCube(:));

    % get dvhPoints for every structure and every scenario the same
    n = 1000;
    if DVHType == 1
        dvhPoints = linspace(0,maxDose*1.05,n);
    else
        dvhPoints = linspace(minDose * 1.05,maxDose*1.05,n);
    end
end

numOfVois = size(cst,1);
dvhV{numOfVois} = [];

for j = 1:numOfVois
    dvhV{j} = getDVHPoints(cst, j, doseCube, dvhPoints, DVHType);
end

% output:
% dvhPoints
% dvhV

end

function [dvh] = getDVHPoints(cst, sIx, doseCube, dvhPoints, DVHType)
n = numel(dvhPoints);
dvh       = NaN * ones(1,n);
% maxDVH = 0;
% counter = 0;
indices     = cst{sIx,4}{1};
numOfVoxels = numel(indices);

doseInVoi   = doseCube(indices);

switch DVHType
    case 1 % cummulative DVH
        for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
        end

    case 2 % differential DVH
        binning = (dvhPoints(2) - dvhPoints(1))/2;
        for j = 1:n % differential DVH        
            dvh(j) = sum(dvhPoints(j) + binning > doseInVoi & doseInVoi > dvhPoints(j) - binning);
        end

end
dvh = dvh ./ numOfVoxels * 100;
end
