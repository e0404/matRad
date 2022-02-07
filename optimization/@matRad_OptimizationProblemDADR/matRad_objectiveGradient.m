function weightGradient = matRad_objectiveGradient(optiProb,wCombined,dij,cst)
% matRad IPOPT callback: gradient function for inverse planning 
% supporting mean dose objectives, EUD objectives, squared overdosage, 
% squared underdosage, squared deviation and DVH objectives
% 
% call
%   g = matRad_gradFuncWrapper(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%
% output
%   g: gradient of objective function
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
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

% get current dose / effect / RBExDose vector
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

% Initializes dose and DADR gradient
doseDADRGradientCombined{1} = zeros(2*dij.doseGrid.numOfVoxels,1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)    
   
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
            
            %Get current optimization function
            objective = cst{i,6}{j};
            
            % only perform gradient computations for objectives
            if isa(objective,'DoseObjectives.matRad_DoseObjective')
                % if we have effect optimization, temporarily replace doses with effect
                if (~isequal(objective.name, 'Mean Dose') && ~isequal(objective.name, 'EUD')) &&...
                    (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection')) 
                    
                    doses = objective.getDoseParameters();
                
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    objective = objective.setDoseParameters(effect);
                end
                
                %dose in VOI
                d_i = d{1}(cst{i,4}{1});
                
                %add to dose gradient
                doseDADRGradientCombined{1}(cst{i,4}{1}) = doseDADRGradientCombined{1}(cst{i,4}{1}) + objective.computeDoseObjectiveGradient(d_i);
            
            % only perform gradient computations for objectives
            elseif isa(objective,'DADRObjectives.matRad_DADRObjective')
                
                %dose-rate in VOI
                index= dij.doseGrid.numOfVoxels *ones([length(cst{i,4}{1}) 1]);
                dadr_i = DADR{1}(cst{i,4}{1});
                
                %add to dose gradient
                doseDADRGradientCombined{1}(index+cst{i,4}{1}) = doseDADRGradientCombined{1}(index+cst{i,4}{1}) + objective.computeDADRObjectiveGradient(dadr_i);                
            end                          
        end
    end
end
  
%project to weight gradient
optiProb.BP_dadr = optiProb.BP_dadr.computeGradient(dij,doseDADRGradientCombined,wCombined);
g = optiProb.BP_dadr.GetGradient();
weightGradient = g{1};


end
