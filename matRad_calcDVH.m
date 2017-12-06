function dvh = matRad_calcDVH(cst,doseCube,dvhType,doseGrid)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad dvh calculation
% 
% call
%   dvh = matRad_calcDVH(cst,doseCube,dvhType,doseGrid)
%
% input
%   cst:                  matRad cst struct
%   doseCube:             arbitrary doseCube (e.g. physicalDose)
%   dvhType: (optional)   string, 'cum' for cumulative, 'diff' for differential
%                         dvh
%   doseGrid: (optional): use predefined evaluation points. Useful when
%                         comparing multiple realizations
%
% output
%   dose volume histogram
%
% References
%   van't Riet et. al., IJROBP, 1997 Feb 1;37(3):731-6.
%   Kataria et. al., J Med Phys. 2012 Oct-Dec; 37(4): 207�213.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('dvhType','var') || isempty(dvhType)
    dvhType = 'cum';
end

if ~exist('doseGrid', 'var') || isempty(doseGrid)
    maxDose = max(doseCube(:));
    minDose = min(doseCube(:));

    % get dvhPoints for every structure and every scenario the same
    n = 1000;
    if strcmp(dvhType, 'cum')
        doseGrid = linspace(0,maxDose*1.05,n);
    elseif strcmp(dvhType, 'diff')
        doseGrid = linspace(0.95*minDose,maxDose*1.05,n);
    end
end

numOfVois = size(cst,1);
dvh = struct;
for i = 1:numOfVois
    dvh(i).doseGrid     = doseGrid;
    dvh(i).volumePoints = getDVHPoints(cst, i, doseCube, doseGrid, dvhType);
    dvh(i).VOIname      = cst{i,2};
end

end %eof 

function dvh = getDVHPoints(cst, sIx, doseCube, dvhPoints, dvhType)
n = numel(dvhPoints);
dvh         = NaN * ones(1,n);
indices     = cst{sIx,4}{1};
numOfVoxels = numel(indices);

doseInVoi   = doseCube(indices);

switch dvhType
    case 'cum' % cummulative DVH
        for j = 1:n
            dvh(j) = sum(doseInVoi >= dvhPoints(j));
        end

    case 'diff' % differential DVH
        binning = (dvhPoints(2) - dvhPoints(1))/2;
        for j = 1:n % differential DVH        
            dvh(j) = sum(dvhPoints(j) + binning > doseInVoi & doseInVoi > dvhPoints(j) - binning);
        end

end
dvh = dvh ./ numOfVoxels * 100;
end %eof getDVHPoints
