function [dvh, qi] = ...
    matRad_calcIndicators(cst,pln,doseCube,dvhType,doseGrid,refGy,refVol,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad indictor calculation (contains QI as well as DVH)
% 
% call
%   matRad_calcIndicators(cst,pln,cube,dvhType,param,refGy,refVol,lineStyleIndicator)
%
% input
%   cst:                  matRad cst struct
%   pln:                  matRad pln struct
%   doseCube:             arbitrary doseCube where (e.g. physicalDose)
%   dvhType: (optional)   string, 'cum' for cumulative, 'diff' for differential
%                         dvh
%   doseGrid: (optional): use predefined evaluation points. Useful when
%                         comparing multiple realizations
%   refGy: (optional)     array of dose values used for V_XGy calculation
%                         default is [40 50 60]
%   refVol:(optional)     array of volumes (0-100) used for D_X calculation
%                         default is [2 5 95 98]
%   param:                parameter struct, e.g. for changing loglevel
%
% output
%   various quality indicators as well as dvh
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
        doseGrid = linspace(minDose * 1.05,maxDose*1.05,n);
    end
end

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end

if ~exist('refVol', 'var') || isempty(refVol)
    refVol = [2 5 98 95];
end

if ~exist('refGy', 'var') || isempty(refGy)
    refGy = linspace(0,max(doseCube(:)),6);
end

dvh = calcDVH(cst, doseCube, dvhType, doseGrid);
qi = matRad_calcQualityIndicators(cst,pln,doseCube,param,refGy,refVol);

    function dvhV = calcDVH(cst, doseCube, dvhType, doseGrid)

    numOfVois = size(cst,1);
    dvhV{numOfVois,1} = [];
    for i = 1:numOfVois
        dvhV{i,1} = [doseGrid; getDVHPoints(cst, i, doseCube, doseGrid, dvhType)];
    end

        function [dvh] = getDVHPoints(cst, sIx, doseCube, dvhPoints, dvhType)
        n = numel(dvhPoints);
        dvh       = NaN * ones(1,n);
        indices     = cst{sIx,4}{1};
        numOfVoxels = numel(indices);

        doseInVoi   = doseCube(indices);

        switch dvhType
            case 'cum' % cummulative DVH
                for j = 1:n
                    dvh(j) = sum(doseInVoi > dvhPoints(j));
                end

            case 'diff' % differential DVH
                binning = (dvhPoints(2) - dvhPoints(1))/2;
                for j = 1:n % differential DVH        
                    dvh(j) = sum(dvhPoints(j) + binning > doseInVoi & doseInVoi > dvhPoints(j) - binning);
                end

        end
        dvh = dvh ./ numOfVoxels * 100;
      end %eof getDVHPoints
    end %eof dvh
end %eof 
