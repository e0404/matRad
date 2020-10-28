function resultGUI = matRad_calcCubes(w,dij,scenNum)
% matRad computation of all cubes for the resultGUI struct 
% which is used as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij)
%   resultGUI = matRad_calcCubes(w,dij,scenNum)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
%
% References
%   -
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

if nargin < 3
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
%

% compute physical dose for all beams individually and together
for i = 1:length(beamInfo)
    resultGUI.(['physicalDose', beamInfo(i).suffix]) = reshape(full(dij.physicalDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
end

% consider RBE for protons
if isfield(dij,'RBE')
   fprintf(['matRad: applying a constant RBE of ' num2str(dij.RBE) ' \n']);
   for i = 1:length(beamInfo)
        resultGUI.(['RBExDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
   end
end

% consider LET
if isfield(dij,'mLETDose')
    for i = 1:length(beamInfo)
        LETDoseCube                                 = dij.mLETDose{scenNum} * (resultGUI.w .* beamInfo(i).logIx);
        resultGUI.(['LET', beamInfo(i).suffix])     = zeros(dij.doseGrid.dimensions);
        ix                                          = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    end
end

if isfield(dij,'physicalDose_MCvar')
    resultGUI.physicalDose_MCvar = reshape(full(dij.physicalDose_MCvar{scenNum} * (resultGUI.w .* beamInfo(i).logIx)),dij.doseGrid.dimensions);
    resultGUI.physicalDose_MCstd = sqrt(resultGUI.physicalDose_MCvar);
    resultGUI.physicalDose_MCstdRel = resultGUI.physicalDose_MCstd ./ resultGUI.physicalDose;
end

% consider biological optimization for carbon ions
if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
   
    ix = dij.bx~=0;

    for i = 1:length(beamInfo)  
       wBeam = (resultGUI.w .* beamInfo(i).logIx);
       resultGUI.(['effect', beamInfo(i).suffix])       = full(dij.mAlphaDose{scenNum} * wBeam + (dij.mSqrtBetaDose{scenNum} * wBeam).^2);
       resultGUI.(['effect', beamInfo(i).suffix])       = reshape(resultGUI.(['effect', beamInfo(i).suffix]),dij.doseGrid.dimensions);
    
       resultGUI.(['RBExDose', beamInfo(i).suffix])     = zeros(size(resultGUI.(['effect', beamInfo(i).suffix])));
       resultGUI.(['RBExDose', beamInfo(i).suffix])(ix) = (sqrt(dij.ax(ix).^2 + 4 .* dij.bx(ix) .* resultGUI.(['effect', beamInfo(i).suffix])(ix)) - dij.ax(ix))./(2.*dij.bx(ix));

       resultGUI.(['RBE', beamInfo(i).suffix])          = resultGUI.(['RBExDose', beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);

       resultGUI.(['alpha', beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);
       resultGUI.(['beta',  beamInfo(i).suffix])        = zeros(dij.doseGrid.dimensions);

       AlphaDoseCube                                    = full(dij.mAlphaDose{scenNum} * wBeam);
       resultGUI.(['alpha', beamInfo(i).suffix])(ix)    = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);

       SqrtBetaDoseCube                                 = full(dij.mSqrtBetaDose{scenNum} * wBeam);
       resultGUI.(['beta', beamInfo(i).suffix])(ix)     = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
    end
end

% group similar fields together
resultGUI = orderfields(resultGUI);

% interpolation if dose grid does not match ct grid
if any(dij.ctGrid.dimensions~=dij.doseGrid.dimensions)
   myFields = fieldnames(resultGUI);
   for i = 1:numel(myFields)
      
       if numel(resultGUI.(myFields{i})) == dij.doseGrid.numOfVoxels
           
           % interpolate!
           resultGUI.(myFields{i}) = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                             resultGUI.(myFields{i}), ...
                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);
           
       end
       
   end   
end

end



