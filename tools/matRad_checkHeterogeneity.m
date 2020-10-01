function pln = matRad_checkHeterogeneity(pln,cst,baseData)
% script to check the consitency of heterogeneity parameters and setting 
% default parameters for heterogeneity correction
%
% call
%   pln = matRad_checkHeterogeneity(pln,cst)
%
% input
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   baseData:       matRad machine base data
%
% output
%   pln:            updated matRad plan meta information struct
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

if isfield(pln,'heterogeneity')
    if ~isfield(pln.heterogeneity,'useDoseCurves') || isempty(pln.heterogeneity.useDoseCurves)
        pln.heterogeneity.useDoseCurves = false;
    end
    if ~isfield(pln.heterogeneity,'calcHetero') || isempty(pln.heterogeneity.calcHetero)
        pln.heterogeneity.calcHetero = false;
    end
    if ~isfield(pln.heterogeneity,'type') || isempty(pln.heterogeneity.type)
        if pln.heterogeneity.calcHetero
            pln.heterogeneity.type = 'complete';
            warning('Heterogeneity correction wanted but no correction type specified, default was selected.'); 
        else
            pln.heterogeneity.type = 'none';
        end
    end
else
     pln.heterogeneity.useDoseCurves = false;
     pln.heterogeneity.calcHetero = false; 
     pln.heterogeneity.type = 'none';
end

if pln.heterogeneity.calcHetero
    fprintf('Heterogeneity correction enabled. \n');
    heteroCST = false;
    for i = 1:length(cst(:,1))
        if isfield(cst{i,5},'HeterogeneityCorrection')
            heteroCST = true;
            break
        end
    end
    if ~isstruct(baseData.Z) || ~heteroCST
       warning('Heterogeneity correction enabled but no usable data in cst or unsuitable base data. Correction cannot be applied.'); 
       pln.heterogeneity.calcHetero = false;
    end
else
    fprintf('Heterogeneity correction disabled. \n');
end

end