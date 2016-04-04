function resultGUI = matRad_calcCubes(w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad computation of all cubes for the resultGUI struct which is used
% as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij,cst)
%
% input
%   w:   bixel weight vector
%   dij: dose influence matrix
%   cst: matRad cst struct
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

resultGUI.w = w;

% calc dose and reshape from 1D vector to 2D array
resultGUI.physicalDose = reshape(dij.physicalDose*resultGUI.w,dij.dimensions);

if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')

    a_x = zeros(size(resultGUI.physicalDose));
    b_x = zeros(size(resultGUI.physicalDose));

    for  i = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') 
            a_x(cst{i,4}) = cst{i,5}.alphaX;
            b_x(cst{i,4}) = cst{i,5}.betaX;
        end
    end
    
    resultGUI.effect = (dij.mAlphaDose*resultGUI.w+(dij.mSqrtBetaDose*resultGUI.w).^2);
    resultGUI.effect = reshape(resultGUI.effect,dij.dimensions);
    
    resultGUI.RBExDose     = zeros(size(resultGUI.effect));
    ix                     = resultGUI.effect>0;
    resultGUI.RBExDose(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.effect(ix)) - a_x(ix))./(2.*b_x(ix)));
    resultGUI.RBE          = resultGUI.RBExDose./resultGUI.physicalDose;
   
    AlphaDoseCube    = dij.mAlphaDose * resultGUI.w;
    resultGUI.alpha  = (reshape(AlphaDoseCube,dij.dimensions))./resultGUI.physicalDose;
    SqrtBetaDoseCube = dij.mSqrtBetaDose * resultGUI.w;
    resultGUI.beta   = ((reshape(SqrtBetaDoseCube,dij.dimensions))./resultGUI.physicalDose).^2;
    
end