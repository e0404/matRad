function resultGUI = matRad_calcCubes(w,dij,cst,scenNum,pln)
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
resultGUI.physicalDose = reshape(dij.physicalDose{scenNum}*resultGUI.w,dij.dimensions);

% consider RBE for protons
if isfield(dij,'RBE')
   fprintf(['matRad: applying a constant RBE of ' num2str(dij.RBE) ' \n']); 
   resultGUI.RBExDose     = resultGUI.physicalDose * dij.RBE;
end

% consider VOI priorities
[cst,resultGUI.overlapCube]  = matRad_setOverlapPriorities(cst,dij.dimensions);

if isfield(dij,'mLETDose')
    LETDoseCube       = dij.mLETDose{scenNum} * resultGUI.w;
    resultGUI.LET     = zeros(dij.dimensions);
    ix                = resultGUI.physicalDose>0;
    resultGUI.LET(ix) = LETDoseCube(ix)./resultGUI.physicalDose(ix);
end

% consider biological optimization for carbon ions
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


% if linear scaling model was used - check if the used lambda makes sense
if isfield(pln.bioParam,'lamda_1_1')
   
   % find all target voxels from cst cell array
   targetAlphaX = [];
   targetBetaX  = [];
   numVoxel     = [];
   V = [];
   for i=1:size(cst,1)
       if isequal(cst{i,3},'TARGET') && ( ~isempty(cst{i,6}) || ~isempty(findstr(cst{i,2},'ScenUnion')) )
          targetAlphaX = [targetAlphaX cst{i,5}.alphaX]; targetBetaX  = [targetBetaX cst{i,5}.betaX];
          numVoxel     = [numVoxel numel(cst{i,4}{:})];
          V = [V;vertcat(cst{i,4}{:})];
       end
   end

   % normalize alphaX in case of different target tissues
   targetAlphaX = (targetAlphaX .* numVoxel)./(sum(numVoxel));
   targetBetaX  = (targetBetaX .* numVoxel)./(sum(numVoxel));
   % Remove double voxels
   V = unique(V);

   refConstRBE = 1.1;
   meanLET  = mean(resultGUI.LET(V)); 
   meanDose = mean(resultGUI.physicalDose(V));
   meanRBE  = mean(resultGUI.RBE(V));
   
   disp(['mean RBE in target is: ' num2str(meanRBE)]);
   
   meanLambda_1_1 = ( targetAlphaX*(refConstRBE-1) + targetBetaX .* meanDose .*(refConstRBE^2-1) ) ./ (meanLET - pln.bioParam.corrFacEntranceRBE);
   
   relDiff =  (((meanLambda_1_1/pln.bioParam.lamda_1_1) -1 ) * 100);
   if abs(relDiff) > 10
      warning(['relatve deviation of lamda_1_1: ' num2str(relDiff) '%']);
      warning(['used lamda_1_1: ' num2str(pln.bioParam.lamda_1_1) '  lamda_1_1 of ' num2str(meanLambda_1_1) ' would ensure a mean RBE of 1.1 in the target']);
   end 
      
      %%ToDO: asign non specific normal tissue a RBE of 1.1
end

%% ToDo: calculate the worst case dose for COWC

% write worst case dose distribution if VWWC opt used for one objective
% only implemented for physical dose
saveWCCube = 0;
for i = 1:size(cst,1)
   if ~isempty(cst{i,6})
        for j = 1:numel(cst{i,6})
            if(strcmp(cst{i,6}(j).robustness,'VWWC'))
                saveWCCube = 1;
            end
        end
   end
end

if(saveWCCube == 1)
    quantity = 'physicalDose';
    if isfield(dij,'RBE')
        quantity = 'RBExD';
        for i = 1: length(dij.indexforOpt)
           d{i} = (dij.physicalDose{dij.indexforOpt(i)} * w ) * dij.RBE;
        end
    elseif isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
        quantity = 'RBExD';
        for i = 1: length(dij.indexforOpt)
           effect   = full(dij.mAlphaDose{dij.indexforOpt(i)}*resultGUI.w+(dij.mSqrtBetaDose{dij.indexforOpt(i)}*resultGUI.w).^2);
           d{i}     = zeros(size(effect));
           ix       = effect>0;
           d{i}(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* effect(ix)) - a_x(ix))./(2.*b_x(ix)));
        end
    else
        for i = 1: length(dij.indexforOpt)
           d{i} = dij.physicalDose{dij.indexforOpt(i)} * w;
        end
    end
    
    [d_max,~] = max([d{:}],[],2);  
    [d_min,~] = min([d{:}],[],2);
             
    d_wc = d_max;
    for  i = 1:size(cst,1)
        for j = 1:numel(cst{i,6})
  
            if isequal(cst{i,3},'TARGET')
                d_wc(cst{i,4}{1}) = d_min(cst{i,4}{1});
            end
        end
    end
    resultGUI.(['VWWC_' quantity]) = reshape(d_wc,dij.dimensions);
end
