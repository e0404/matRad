function f = matRad_objectiveFunction(optiProb,wCombined,dij,cst)
% matRad IPOPT objective function wrapper
% 
% call
%   f = matRad_objectiveFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: matRad optimization problem
%   w:        beamlet/ pencil beam weight vector
%   dij:      matRad dose influence struct
%   cst:      matRad cst struct
%
% output
%   f: objective function value
%
% References
%   -
%
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

% get current dose / effect / RBExDose vector
%d = optiProb.matRad_backProjection(w,dij);

indexPB = numel(wCombined)/2;
w = wCombined(1:indexPB);
I = wCombined(indexPB+1:end);

%Compute dose rate
optiProb.BP_dadr = optiProb.BP_dadr.compute(dij,wCombined);
dCombined = optiProb.BP_dadr.GetResult();

optiProb.BP = optiProb.BP.compute(dij,w);
dLegacy = optiProb.BP.GetResult();

%separate into dose & dose-rate
index = dij.doseGrid.numOfVoxels;
d{1} = dCombined{1}(1:index);
DADR{1} =dCombined{1}(index+1:end);


% Initialize f
f = 0;

% compute objectiveective function for every VOI.
for  i = 1:size(cst,1)    
    if ~isempty(cst{i,4}{1}) &&( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
              objective = cst{i,6}{j};
             if isa(objective,'DoseObjectives.matRad_DoseObjective') 
                 d_i = d{1}(cst{i,4}{1});                 
                 f = f + objective.computeDoseObjectiveFunction(d_i);
             elseif isa(objective,'DADRObjectives.matRad_DADRObjective')
                 dadr_i = DADR{1}(cst{i,4}{1});                 
                 f = f + objective.computeDADRObjectiveFunction(dadr_i);
             end
        end                          
    end  
end
