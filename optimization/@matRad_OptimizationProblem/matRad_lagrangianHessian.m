function hess = matRad_lagrangianHessian(optiProb,w,dij,cst,sigma,lambda)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: hessian of the lagrangian function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint
%
% call
%   jacob = matRad_lagrangianHessian(w,dij,cst,optiProb,sigma,lambda)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   lambda: lagrangian multipliers for constraints
%   sigma:  multipliers for objectives
%
% output
%   hess: hessian of the lagrangian
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

hessDose = sparse(numel(d{1}),numel(d{1}));

lambdaCounter = 1;

% compute objective function for every VOI.
for i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        d_i = d{1}(cst{i,4}{1});
        
        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            obj = cst{i,6}{j}; %Get the Optimization Object
            
            %Default empty hessian
            hessianTmp = sparse(numel(d_i),numel(d_i));
            
            % only perform computations for constraints
            %if ~isempty(strfind(obj.type,'constraint'))
            if isa(obj,'DoseConstraints.matRad_DoseConstraint')
                
                % compute reference
                if isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection')
                    
                    doses = obj.getDoseParameters();
                    
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    obj = obj.setDoseParameters(effect);
                end
                
                %Select the subset of lambdas that 
                nCfuncs = obj.numConstraints(numel(d_i));
                lambdaSub = lambda(lambdaCounter:lambdaCounter+nCfuncs-1);
                                
                hessianTmp = obj.computeDoseConstraintHessian(lambdaSub,d_i);
                
                lambdaCounter = lambdaCounter + 1;
            end
            
            if isa(obj,'DoseObjectives.matRad_DoseObjective')
                % compute reference
                if isa(optiProb.BP,'matRad_EffectProjection') && ~isa(optiProb.BP,'matRad_VariableRBEProjection')
                    
                    doses = obj.getDoseParameters();
                    
                    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                    obj = obj.setDoseParameters(effect);
                end

                hessianTmp = sigma*obj.computeDoseObjectiveHessian(d_i);
            end
            
            hessDose(cst{i,4}{1},cst{i,4}{1}) = hessDose(cst{i,4}{1},cst{i,4}{1}) + hessianTmp;                                                
        end
        
    end
    
end

%Hessian construction
for i =1:dij.numOfScenarios    
    if isa(optiProb.BP,'matRad_DoseProjection') || isa(optiPorb.BP,'matRad_ConstantRBEProjection')
        hess = sparse(tril(dij.physicalDose{i}'*hessDose*dij.physicalDose{i}));        
    elseif isa(optiProb.BP,'matRad_EffectProjection')
        hess = sparse(tril(dij.mSqrtBetaDose{i}'*hessDose*dij.physicalDose{i}));        
        
    else
        error('Unknown Dose Projection!');
    end   
end

%{
% Calculate jacobian with dij projections
for i = 1:dij.numOfScenarios

    elseif isa(optiProb.BP,'matRad_EffectProjection')
        
        if ~isempty(mSqrtBetaDoseProjection) && ~isempty(mAlphaDoseProjection)
            
            jacobLogical            = (scenID == i);
            jacobLogical2           = (scenID2 == i);
            mSqrtBetaDoseProjection = mSqrtBetaDoseProjection(:,jacobLogical2)' * dij.mSqrtBetaDose{i} * w;
            mSqrtBetaDoseProjection = sparse(voxelID(jacobLogical2),constraintID(jacobLogical2),mSqrtBetaDoseProjection,...
                size(mAlphaDoseProjection(:,jacobLogical),1),size(mAlphaDoseProjection(:,jacobLogical),2));
            
            hess(jacobLogical,:)   = mAlphaDoseProjection(:,jacobLogical)' * dij.mAlphaDose{i} +...
                mSqrtBetaDoseProjection' * dij.mSqrtBetaDose{i};
            
        end
    end
end
%}

end
