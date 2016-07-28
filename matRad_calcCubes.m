function resultGUI = matRad_calcCubes(w,dij,cst,scenNum)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad computation of all cubes for the resultGUI struct which is used
% as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij,cst)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
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

if nargin < 4
    scenNum = 1;
end

resultGUI.w = w;

% calc dose and reshape from 1D vector to 2D array
resultGUI.physicalDose = reshape(full(dij.physicalDose{scenNum}*resultGUI.w),dij.dimensions);

% consider VOI priorities
[cst,resultGUI.overlapCube]  = matRad_setOverlapPriorities(cst,dij.dimensions);

if isfield(dij,'mLETDose')
    LETDoseCube       = dij.mLETDose{scenNum} * resultGUI.w;
    resultGUI.LET     = zeros(dij.dimensions);
    ix                = resultGUI.physicalDose>0;
    resultGUI.LET(ix) = LETDoseCube(ix)./resultGUI.physicalDose(ix);

end

if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')

    a_x = zeros(size(resultGUI.physicalDose));
    b_x = zeros(size(resultGUI.physicalDose));

    for i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}{scenNum}) = cst{i,5}.alphaX;
            b_x(cst{i,4}{scenNum}) = cst{i,5}.betaX;
        end
    end
    
    resultGUI.effect = full(dij.mAlphaDose{scenNum}*resultGUI.w+(dij.mSqrtBetaDose{scenNum}*resultGUI.w).^2);
    resultGUI.effect = reshape(resultGUI.effect,dij.dimensions);
    
    resultGUI.RBExDose     = zeros(size(resultGUI.effect));
    ix                     = resultGUI.effect>0;
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));
    resultGUI.RBE          = resultGUI.RBExDose./resultGUI.physicalDose;
   
    resultGUI.alpha     = zeros(size(resultGUI.effect));
    resultGUI.beta      = zeros(size(resultGUI.effect));
    AlphaDoseCube       = full(dij.mAlphaDose{scenNum} * resultGUI.w);
    resultGUI.alpha(ix) = AlphaDoseCube(ix)./resultGUI.physicalDose(ix);
    SqrtBetaDoseCube    = full(dij.mSqrtBetaDose{scenNum} * resultGUI.w);
    resultGUI.beta(ix)  = (SqrtBetaDoseCube(ix)./resultGUI.physicalDose(ix)).^2;
    
end
