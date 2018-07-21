function resultGUI = matRad_calcCubes(w,dij,cst,scenNum,calcBeamDose)
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
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('scenNum','var') || isempty(scenNum)
    scenNum = 1;
end

if ~exist('calcBeamDose','var') || isempty(calcBeamDose)
    calcBeamDose = true;
end

resultGUI.w = w;

% indicating the nominal and expected dose influence matrix
if scenNum > 1
    probQuant = {''};
else
    probQuant = {'','Exp'};   
end

% get bixel - beam correspondence 
if calcBeamDose
    for i = 1:dij.numOfBeams
        beamInfo(i).suffix = ['_beam', num2str(i)];
        beamInfo(i).logIx  = (dij.beamNum == i);      
    end
    beamInfo(dij.numOfBeams+1).suffix = '';
    beamInfo(dij.numOfBeams+1).logIx  = true(size(w));
else
    beamInfo(1).suffix = '';
    beamInfo(1).logIx  = true(size(w));
end

% compute physical dose for all beams individually and together
for j = 1:length(probQuant)
    for i = 1:length(beamInfo)
        resultGUI.(['physicalDose' probQuant{1,j} beamInfo(i).suffix]) = reshape(full(dij.(['physicalDose' probQuant{1,j}]){scenNum} *...
                                                                           (resultGUI.w .* beamInfo(i).logIx)),dij.dimensions);                                                            
    end
end

% consider RBE for protons
for j = 1:length(probQuant)
    if isfield(dij,'RBE')
       for i = 1:length(beamInfo)
            resultGUI.(['RBExD' probQuant{1,j} beamInfo(i).suffix]) = resultGUI.(['physicalDose' probQuant{1,j} beamInfo(i).suffix]) * dij.RBE;
       end
    end
end

% consider VOI priorities
[cst,resultGUI.overlapCube]  = matRad_setOverlapPriorities(cst,dij.dimensions);

% consider LET
for j = 1:length(probQuant)
   if isfield(dij,['mLETDose' probQuant{1,j}]) && scenNum == 1
       for i = 1:length(beamInfo)
           LETDoseCube                                 = dij.(['mLETDose' probQuant{1,j}]){scenNum} * (resultGUI.w .* beamInfo(i).logIx);
           resultGUI.(['LET' probQuant{1,j} beamInfo(i).suffix])     = zeros(dij.dimensions);
           ix                                          = resultGUI.(['physicalDose' probQuant{1,j} beamInfo(i).suffix]) > 0;
           resultGUI.(['LET' probQuant{1,j} beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose' probQuant{1,j} beamInfo(i).suffix])(ix);
       end
   end
end

for j = 1:length(probQuant)
   % consider biological optimization for carbon ions
   if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')

       a_x = zeros(dij.dimensions);
       b_x = zeros(dij.dimensions);

       for l = 1:size(cst,1)
           % Only take OAR or target VOI.
           if isequal(cst{l,3},'OAR') || isequal(cst{l,3},'TARGET') 
               a_x(cst{l,4}{1}) = cst{l,5}.alphaX;
               b_x(cst{l,4}{1}) = cst{l,5}.betaX;
           end
       end

       ix = b_x~=0;

       for i = 1:length(beamInfo)  
          wBeam = (resultGUI.w .* beamInfo(i).logIx);
          resultGUI.(['effect' probQuant{1,j} beamInfo(i).suffix])     = full(dij.(['mAlphaDose' probQuant{1,j}]){scenNum} * wBeam ...
                                                                           + (dij.(['mSqrtBetaDose' probQuant{1,j}]){scenNum} * wBeam).^2);
          resultGUI.(['effect' probQuant{1,j} beamInfo(i).suffix])     = reshape(resultGUI.(['effect' probQuant{1,j} beamInfo(i).suffix]),dij.dimensions);

          resultGUI.(['RBExD' probQuant{1,j} beamInfo(i).suffix])      = zeros(size(resultGUI.(['effect' probQuant{1,j} beamInfo(i).suffix])));
          resultGUI.(['RBExD' probQuant{1,j} beamInfo(i).suffix])(ix)  = (sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.(['effect' probQuant{1,j} beamInfo(i).suffix])(ix)) - a_x(ix))./(2.*b_x(ix));

          resultGUI.(['RBE' probQuant{1,j} beamInfo(i).suffix])        = resultGUI.(['RBExD' probQuant{1,j} beamInfo(i).suffix])./resultGUI.(['physicalDose' probQuant{1,j} beamInfo(i).suffix]);

          resultGUI.(['alpha', beamInfo(i).suffix])        = zeros(dij.dimensions);
          resultGUI.(['beta',  beamInfo(i).suffix])        = zeros(dij.dimensions);

          AlphaDoseCube                                    = full(dij.mAlphaDose{scenNum} * wBeam);
          resultGUI.(['alpha', beamInfo(i).suffix])(ix)    = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);

          SqrtBetaDoseCube                                 = full(dij.mSqrtBetaDose{scenNum} * wBeam);
          resultGUI.(['beta', beamInfo(i).suffix])(ix)     = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
       end
   end
end
% group similar fields together
resultGUI = orderfields(resultGUI);

% rename fields
robID = '';
for i = 1:size(cst,1)  
    for j = 1:numel(cst{i,6})
        if ~strcmp(cst{i,6}(j).robustness,'none')
            robID = 'Rob'; break;
        end
    end   
end

fields = fieldnames(resultGUI);
for i=1:numel(fields)
    if ~strcmp(fields{i,1},'overlapCube') && ~isempty(robID)
        resultGUI.([fields{i,1} robID]) =  resultGUI.(fields{i,1});
        resultGUI = rmfield(resultGUI,fields{i,1});
    end 
end


end



