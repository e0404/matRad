function resultGUI = matRad_calcCubes(w,dij,cst,scenNum)
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

% get bixel - beam correspondence  
for i = 1:dij.numOfBeams
    beamInfo(i).suffix = ['_beam', num2str(i)];
    beamInfo(i).logIx  = (dij.beamNum == i);      
end
beamInfo(dij.numOfBeams+1).suffix = '';
beamInfo(dij.numOfBeams+1).logIx  = true(size(w));

% compute physical dose for all beams individually and together
for i = 1:length(beamInfo)
    resultGUI.(['physicalDose', beamInfo(i).suffix]) = reshape(full(dij.physicalDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.dimensions);
end

% consider RBE for protons
if isfield(dij,'RBE')
   for i = 1:length(beamInfo)
        resultGUI.(['RBExD', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
   end
end

% consider VOI priorities
[cst,resultGUI.overlapCube]  = matRad_setOverlapPriorities(cst,dij.dimensions);

% consider LET
if isfield(dij,'mLETDose')
    for i = 1:length(beamInfo)
        LETDoseCube                                 = dij.mLETDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx);
        resultGUI.(['LET', beamInfo(i).suffix])     = zeros(dij.dimensions);
        ix                                          = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    end
end


% consider biological optimization for carbon ions
if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
   
    for i = 1:length(beamInfo)  
   
       wBeam = (resultGUI.w .* beamInfo(i).logIx);
       
       ix = dij.betaX~=0 & resultGUI.(['physicalDose', beamInfo(i).suffix])(:) > 0;

       resultGUI.(['effect', beamInfo(i).suffix])       = full(dij.mAlphaDose{scenNum} * wBeam + (dij.mSqrtBetaDose{scenNum} * wBeam).^2);
       resultGUI.(['effect', beamInfo(i).suffix])       = reshape(resultGUI.(['effect', beamInfo(i).suffix]),dij.dimensions);
    
       resultGUI.(['RBExD', beamInfo(i).suffix])        = zeros(size(resultGUI.(['effect', beamInfo(i).suffix])));
       resultGUI.(['RBExD', beamInfo(i).suffix])(ix)    = (sqrt(dij.alphaX(ix).^2 + 4 .* dij.betaX(ix) .* resultGUI.(['effect', beamInfo(i).suffix])(ix)) - dij.alphaX(ix))./(2.*dij.betaX(ix));

       resultGUI.(['RBE', beamInfo(i).suffix])          = resultGUI.(['RBExD', beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);

       resultGUI.(['alpha', beamInfo(i).suffix])        = zeros(dij.dimensions);
       resultGUI.(['beta',  beamInfo(i).suffix])        = zeros(dij.dimensions);

       AlphaDoseCube                                    = full(dij.mAlphaDose{scenNum} * wBeam);
       resultGUI.(['alpha', beamInfo(i).suffix])(ix)    = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);

       SqrtBetaDoseCube                                 = full(dij.mSqrtBetaDose{scenNum} * wBeam);
       resultGUI.(['beta', beamInfo(i).suffix])(ix)     = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
    end
end

% group similar fields together
resultGUI = orderfields(resultGUI);

end



