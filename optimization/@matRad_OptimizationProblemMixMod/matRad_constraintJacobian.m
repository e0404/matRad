function jacob = matRad_constraintJacobian(optiProb,w,dij,cst)
% matRad IPOPT callback: jacobian function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%
% output
%   jacob: jacobian of constraint function
%
% References
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
%d = matRad_backProjection(w,dij,optiProb);
%d = optiProb.matRad_backProjection(w,dij);
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% initialize jacobian (only single scenario supported in optimization)
jacob = sparse([]);

% initialize projection matrices and id containers
DoseProjection{1}          = sparse([]);
mAlphaDoseProjection{1}    = sparse([]);
mSqrtBetaDoseProjection{1} = sparse([]);
voxelID                     = [];
constraintID                = [];

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
fullScen      = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
contourScen   = fullScen{1};

% compute objective function for every VOI.
for i = 1:size(cst,1)
   
   % Only take OAR or target VOI.
   if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
      
      % loop over the number of constraints for the current VOI
      for j = 1:numel(cst{i,6})
         
         constraint = cst{i,6}{j}; %Get the Optimization Object
         
         % only perform computations for constraints
         if isa(constraint,'DoseConstraints.matRad_DoseConstraint')
            
            % retrieve the robustness type
            robustness = constraint.robustness;
            
            % rescale dose parameters to biological optimization quantity if required
            constraint = optiProb.BP.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            switch robustness
               
               case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                     d_i = d{1}(cst{i,4}{1});
                     jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               case 'PROB' % if prob opt: sum up expectation value of objectives
                  
                  d_i = dExp{1}(cst{i,4}{1});
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~exist('d_tmp','var')
                     d_tmp = [d{useScen}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_max;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_min;
                  end
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               case 'VWWC_INV'  %inverse voxel-wise conformitiy - takes maximum dose in TARGET and minimum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~exist('d_tmp','var')
                     d_tmp = [d{useScen}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_min;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_max;
                  end
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
            end
            
            nConst = size(jacobSub,2);
            
            %Iterate through columns of the sub-jacobian
            if isa(optiProb.BP,'matRad_DoseProjection') && ~isempty(jacobSub) || isa(optiProb.BP,'matRad_ConstantRBEProjection')
               
               startIx = size(DoseProjection{1},2) + 1;
               %First append the Projection matrix with sparse zeros
               DoseProjection{1}          = [DoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];
               
               %Now directly write the jacobian in there
               DoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;
               
            elseif isa(optiProb.BP,'matRad_EffectProjection') && ~isempty(jacobSub)
               
               if isa(optiProb.BP,'matRad_VariableRBEProjection')
                  scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
                  jacobSub     = jacobSub./(2*dij.bx(cst{i,4}{1}) .* scaledEffect);
               end
               
               startIx = size(mAlphaDoseProjection{1},2) + 1;
               
               %First append the alphaDose matrix with sparse
               %zeros then insert
               mAlphaDoseProjection{1}    = [mAlphaDoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];
               mAlphaDoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;
               
               %The betadose has a different structure due to the
               %quadratic transformation, but in principle the
               %same as above
               mSqrtBetaDoseProjection{1} =  [mSqrtBetaDoseProjection{1}, sparse(repmat(cst{i,4}{1},nConst,1),repmat(1:numel(cst{i,4}{1}),1,nConst),2*reshape(jacobSub',[],1),dij.doseGrid.numOfVoxels,nConst*numel(cst{i,4}{1}))];
               
               if isempty(constraintID)
                  newID = 1;
               else
                  newID = constraintID(end)+1;
               end
               
               voxelID = [voxelID;repmat(cst{i,4}{1},nConst,1)];                         %Keep track of voxels for organizing the sqrt(beta)Dose projection later
               constraintID = [constraintID, ...
                  reshape(ones(numel(cst{i,4}{1}),1)*[newID:newID+nConst-1],[1 nConst*numel(cst{i,4}{1})])];  %Keep track of constraints for organizing the sqrt(beta)Dose projection later
            end
            
        
         end
         
      end
      
   end
   
end


scenario = 1;
% enter if statement also for protons using a constant RBE
if isa(optiProb.BP,'matRad_DoseProjection')
   
   if ~isempty(DoseProjection{scenario})
      jacob = DoseProjection{scenario}' * dij.physicalDose{scenario};
   end
   
elseif isa(optiProb.BP,'matRad_ConstantRBEProjection')
   
   if ~isempty(DoseProjection{scenario})
      jacob = DoseProjection{scenario}' * dij.RBE * dij.physicalDose{scenario};
   end
   
elseif isa(optiProb.BP,'matRad_EffectProjection')
   
   if ~isempty(mSqrtBetaDoseProjection{scenario}) && ~isempty(mAlphaDoseProjection{scenario})
      mSqrtBetaDoseProjection{scenario} = mSqrtBetaDoseProjection{scenario}' * dij.mSqrtBetaDose{scenario} * w;
      mSqrtBetaDoseProjection{scenario} = sparse(voxelID,constraintID,mSqrtBetaDoseProjection{scenario},...
         size(mAlphaDoseProjection{scenario},1),size(mAlphaDoseProjection{scenario},2));
      
      jacob   = mAlphaDoseProjection{scenario}' * dij.mAlphaDose{scenario} +...
         mSqrtBetaDoseProjection{scenario}' * dij.mSqrtBetaDose{scenario};
      
   end
end

