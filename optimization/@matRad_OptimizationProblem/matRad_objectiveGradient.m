function weightGradient = matRad_objectiveGradient(optiProb,w,dij,cst)
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

matRad_cfg = MatRad_Config.instance();

% get current dose / effect / RBExDose vector
optiProb.BP = optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% also get probabilistic quantities (nearly no overhead if empty)
[dExp,dOmega] = optiProb.BP.GetResultProb();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
fullScen      = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
contourScen   = fullScen{1};

doseGradient          = cell(size(dij.physicalDose));
doseGradient(useScen) = {zeros(dij.doseGrid.numOfVoxels,1)};

%For probabilistic optimization
vOmega = 0;

%For COWC
f_COWC = zeros(size(dij.physicalDose));

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
            
            % retrieve the robustness type
            robustness = objective.robustness;
            
            % rescale dose parameters to biological optimization quantity if required
            objective = optiProb.BP.setBiologicalDosePrescriptions(objective,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            switch robustness
               case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                  for s = 1:numel(useScen)
                     ixScen = useScen(s);
                     ixContour = contourScen(s);
                     d_i = d{ixScen}(cst{i,4}{ixContour});
                     %add to dose gradient
                     doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + objective.computeDoseObjectiveGradient(d_i);
                  end
               case 'STOCH' % perform stochastic optimization with weighted / random scenarios
                  for s = 1:numel(useScen)
                     ixScen = useScen(s);
                     ixContour = contourScen(s);
                     
                     d_i = d{ixScen}(cst{i,4}{ixContour});
                     
                     doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + ...
                        (objective.computeDoseObjectiveGradient(d_i) * scenProb(s));
                     
                  end
                  
               case 'PROB' % use the expectation value and the integral variance influence matrix
                  %First check the speficic cache for probabilistic
                  %opt
                  if ~exist('doseGradientExp','var')
                     doseGradientExp{1} = zeros(dij.doseGrid.numOfVoxels,1);
                  end
                  
                  d_i = dExp{1}(cst{i,4}{1});
                  
                  doseGradientExp{1}(cst{i,4}{1}) = doseGradientExp{1}(cst{i,4}{1}) + objective.computeDoseObjectiveGradient(d_i);
                  
                  p = objective.penalty/numel(cst{i,4}{1});
                  if j == 1
                     vOmega = vOmega + p * dOmega{i,1};
                  end
                  
               case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector for voxel-wise worst case
                  if ~exist('d_tmp','var')
                     d_tmp = [d{useScen}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  [d_max,max_ix] = max(d_Scen,[],2);
                  [d_min,min_ix] = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_max;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_min;
                  end
                  
                  if any(isnan(d_i))
                     matRad_cfg.dispWarning('%d NaN values in gradient.',numel(isnan(d_i)));
                  end
                  
                  deltaTmp = objective.computeDoseObjectiveGradient(d_i);
                  
                  for s = 1:numel(useScen)
                     ixScen = useScen(s);
                     ixContour = contourScen(s);
                     
                     if isequal(cst{i,3},'OAR')
                        currWcIx = double(max_ix == s);
                     elseif isequal(cst{i,3},'TARGET')
                        currWcIx = double(min_ix == s);
                     end
                     
                     doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp.*currWcIx;
                  end
                  
                  
               case 'VWWC_INV'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector for voxel-wise worst case
                  if ~exist('d_tmp','var')
                     d_tmp = [d{:}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  [d_max,max_ix] = max(d_Scen,[],2);
                  [d_min,min_ix] = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_max;
                     matRad_cfg.dispWarning('VWWC_INV is typically a target objective only');
                  elseif isequal(cst{i,3},'TARGET') && isequal(objective.name,'Squared Overdosing')
                     d_i = d_max;
                  elseif isequal(cst{i,3},'TARGET') && isequal(objective.name,'Squared Underdosing')
                     d_i = d_min;
                  elseif isequal(cst{i,3},'TARGET') && isequal(objective.name,'Squared Deviation')  % minimize the minium worst case dose and the maximum worst case dose
                     d_i = [d_max d_min];
                  else
                     matRad_cfg.dispError('Configuration not supported');
                  end
                  
                  if any(isnan(d_i))
                     matRad_cfg.dispWarning('%d NaN values in gradFuncWrapper.',numel(isnan(d_i)));
                  end
                  
                  deltaTmp = objective.computeDoseObjectiveGradient(d_i);
                  
                  for s = 1:numel(useScen)
                     ixScen = useScen(s);
                     ixContour = contourScen(s);
                     
                     if isequal(cst{i,3},'OAR')
                        currWcIx = double(max_ix == s);
                     elseif isequal(cst{i,3},'TARGET')  && isequal(objective.name,'Squared Overdosing')
                        currWcIx = double(max_ix == s);
                     elseif isequal(cst{i,3},'TARGET')  && isequal(objective.name,'Squared Underdosing')
                        currWcIx = double(min_ix == s);
                     elseif isequal(cst{i,3},'TARGET')  && isequal(objective.name,'Squared Deviation')
                        currWcIx = double(max_ix == s) ;
                        currWcIxMin = double(min_ix == s) ;
                     end
                     
                     doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp(:,1).*currWcIx;
                     
                     if isequal(cst{i,3},'TARGET')  && isequal(objective.name,'Squared Deviation')
                        doseGradient{ixScen}(cst{i,4}{ixContour}) = doseGradient{ixScen}(cst{i,4}{ixContour}) + deltaTmp(:,2).*currWcIxMin;
                     end
                  end
                  
               case 'COWC' % composite worst case consideres ovarall the worst objective function value
                  %First check the speficic cache for COWC
                  if ~exist('delta_COWC','var')
                     delta_COWC         = cell(size(doseGradient));
                     delta_COWC(useScen)    = {zeros(dij.doseGrid.numOfVoxels,1)};
                  end
                  
                  for s = 1:numel(useScen)
                     ixScen = useScen(s);
                     ixContour = contourScen(s);
                     
                     d_i = d{ixScen}(cst{i,4}{ixContour});
                     f_COWC(ixScen)                          = f_COWC(ixScen) + objective.computeDoseObjectiveFunction(d_i);
                     delta_COWC{ixScen}(cst{i,4}{ixContour}) = delta_COWC{ixScen}(cst{i,4}{ixContour}) + objective.computeDoseObjectiveGradient(d_i);
                  end
                  
               case 'OWC' % objective-wise worst case consideres the worst individual objective function value
                  %First check the speficic cache for COWC
             
                  f_OWC = zeros(numel(useScen),1);
                  
                  for s = 1:numel(useScen)
                     ixScen    = useScen(s);
                     ixContour = contourScen(s);
                     
                     d_i = d{ixScen}(cst{i,4}{ixContour});
                     f_OWC(s) = objective.computeDoseObjectiveFunction(d_i);
                  end
                  
                  [~,ix] = max(f_OWC);
                  
                   d_i = d{useScen(ix)}(cst{i,4}{contourScen(ix)});
                   
                   doseGradient{1}(cst{i,4}{contourScen(ix)}) = doseGradient{1}(cst{i,4}{contourScen(ix)}) + objective.computeDoseObjectiveGradient(d_i);
                  
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not supported!',objective.robustness);
                  
            end
         end
      end
   end
end

if exist('delta_COWC','var')
   [~,ixCurrWC]           = max(f_COWC(:));
   doseGradient{ixCurrWC} = delta_COWC{ixCurrWC};
end

weightGradient = zeros(dij.totalNumOfBixels,1);

optiProb.BP = optiProb.BP.computeGradient(dij,doseGradient,w);
g = optiProb.BP.GetGradient();

for s = 1:numel(useScen)
   weightGradient = weightGradient + g{useScen(s)};
end

if vOmega ~= 0
   optiProb.BP = optiProb.BP.computeGradientProb(dij,doseGradientExp,vOmega,w);
   gProb = optiProb.BP.GetGradientProb();
   
   %Only implemented for first scenario now
   weightGradient = weightGradient + gProb{1};
end
