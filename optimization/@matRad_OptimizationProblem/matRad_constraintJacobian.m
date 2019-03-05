function jacob = matRad_constraintJacobian(optiProb,w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(w,dij,cst,optiProb)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   optiProb: option struct defining the type of optimization
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
optiProb.BP = optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% initialize jacobian
jacob = sparse([]);

% initialize projection matrices and id containers
DoseProjection          = sparse([]);
mAlphaDoseProjection    = sparse([]);
mSqrtBetaDoseProjection = sparse([]);
voxelID                 = [];
constraintID            = 0;
scenID                  = [];
scenID2                 = [];

% compute objective function for every VOI.
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            obj = cst{i,6}{j}; %Get the Optimization Object

            % only perform computations for constraints
            %if ~isempty(strfind(obj.type,'constraint'))
            if isa(obj,'DoseConstraints.matRad_DoseConstraint')
                
                % compute reference
                if (~isequal(obj.name, 'max dose constraint')      && ~isequal(obj.name, 'min dose constraint')      &&...
                    ~isequal(obj.name, 'max mean dose constraint') && ~isequal(obj.name, 'min mean dose constraint') && ...
                    ~isequal(obj.name, 'min EUD constraint')       && ~isequal(obj.name, 'max EUD constraint'))      && ...
                    (isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection'))
                     
                    doses = obj.getDoseParameters();
                
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    obj = obj.setDoseParameters(effect);
                end
                
                % if conventional opt: just add constraints of nominal dose
                %if strcmp(cst{i,6}{j}.robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    %jacobVec =  matRad_jacobFunc(d_i,cst{i,6}{j},d_ref);
                    jacobSub = obj.computeDoseConstraintJacobian(d_i);
                    
                    scenID  = [scenID;ones(size(jacobSub,2),1)];
                    scenID2 = [scenID2;ones(numel(cst{i,4}{1}),1)];
                    
                    %Iterate through columns of the sub-jacobian
                    %TODO: Maybe this could all be function of the projection
                    %Objects???
                    for c = 1:size(jacobSub,2)
                        jacobVec = jacobSub(:,c);
                        
                        if isa(optiProb.BP,'matRad_DoseProjection') && ~isempty(jacobVec) || isa(optiProb.BP,'matRad_ConstantRBEProjection')
                            
                            DoseProjection          = [DoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.doseGrid.numOfVoxels,1)];
                            
                        elseif isequal(optiProb.BP,'matRad_VariableRBEProjection') && ~isempty(jacobVec)
                            
                            scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
                            
                            delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);
                            
                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.doseGrid.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.doseGrid.numOfVoxels,numel(cst{i,4}{1}))];
                            voxelID                 = [voxelID ;cst{i,4}{1}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];
                           
                        
                        elseif isa(optiProb.BP,'matRad_EffectProjection') && ~isempty(jacobVec)
                            
                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.doseGrid.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.doseGrid.numOfVoxels,numel(cst{i,4}{1}))];
                            voxelID                 = [voxelID ;cst{i,4}{1}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];
                                                    
                        end
                    end
                %end

            end

        end

    end

end

if isa(optiProb.BP,'matRad_EffectProjection')
    constraintID = constraintID(2:end);
end

% Calculate jacobian with dij projections
for i = 1:dij.numOfScenarios
   % enter if statement also for protons using a constant RBE
   if isa(optiProb.BP,'matRad_DoseProjection') || isa(optiProb.BP,'matRad_ConstantRBEProjection')

        if ~isempty(DoseProjection)
            
            jacobLogical          = (scenID == i);
            jacob(jacobLogical,:) = DoseProjection(:,jacobLogical)' * dij.physicalDose{i};
            
        end

    elseif isa(optiProb.BP,'matRad_EffectProjection')

        if ~isempty(mSqrtBetaDoseProjection) && ~isempty(mAlphaDoseProjection)
            
            jacobLogical            = (scenID == i);
            jacobLogical2           = (scenID2 == i);
            mSqrtBetaDoseProjection = mSqrtBetaDoseProjection(:,jacobLogical2)' * dij.mSqrtBetaDose{i} * w;
            mSqrtBetaDoseProjection = sparse(voxelID(jacobLogical2),constraintID(jacobLogical2),mSqrtBetaDoseProjection,...
                                         size(mAlphaDoseProjection(:,jacobLogical),1),size(mAlphaDoseProjection(:,jacobLogical),2));
                                     
            jacob(jacobLogical,:)   = mAlphaDoseProjection(:,jacobLogical)' * dij.mAlphaDose{i} +... 
                                      mSqrtBetaDoseProjection' * dij.mSqrtBetaDose{i};
            
        end
    end
end

end
