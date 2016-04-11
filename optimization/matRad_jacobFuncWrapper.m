function jacob = matRad_jacobFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, exact DVH constraint, max DVH constraint, 
% min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(w,dij,cst,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   jacob: jacobian of constraint function
%
% References
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

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% initialize jacobian
jacob = sparse([]);

% initialize projection matrices and id containers
physicalDoseProjection  = sparse([]);
mAlphaDoseProjection    = sparse([]);
mSqrtBetaDoseProjection = sparse([]);
voxelID                 = [];
constraintID            = 0;

% compute objective function for every VOI.
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if isequal(cst{i,6}(j).type, 'exact DVH constraint') && isequal(type,'effect') || ...
                   isequal(cst{i,6}(j).type, 'max DVH constraint') && isequal(type,'effect') || ...
                   isequal(cst{i,6}(j).type, 'min DVH constraint') && isequal(type,'effect')
                     
                    d_ref = dij.ax(cst{i,4}).*cst{i,6}(j).dose + dij.bx(cst{i,4})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4});

                    jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                    
                    if isequal(type,'none') && ~isempty(jacobVec)

                       physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4},1,jacobVec,dij.numOfVoxels,1)];

                    elseif isequal(type,'effect') && ~isempty(jacobVec)

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4},1,jacobVec,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4},1:numel(cst{i,4}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}))];
                       VoxelID                 = [VoxelID ;cst{i,4}];
                       ConstraintID            = [ConstraintID, repmat(1 + ConstraintID(end),1,numel(cst{i,4}))];

                    elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                       delta = jacobVec./(2*dij.bx(cst{i,4}).*ScaledEffect(cst{i,4}));

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4},1,delta,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4},1:numel(cst{i,4}),2*delta,dij.numOfVoxels,numel(cst{i,4}))];
                       VoxelID                 = [VoxelID ;cst{i,4}];
                       ConstraintID            = [ConstraintID, repmat(1 + ConstraintID(end),1,numel(cst{i,4}))];

                    end

                % if prob opt or voxel-wise worst case: add constraints of all dose scenarios
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
                    
                    for k = 1:dij.numOfScenarios
                        
                        d_i = d{k}(cst{i,4});
                        
                        jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                        
                        
                        if isequal(type,'none') && ~isempty(jacobVec)

                           physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4},1,jacobVec,dij.numOfVoxels,1)];

                        elseif isequal(type,'effect') && ~isempty(jacobVec)

                           mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4},1,jacobVec,dij.numOfVoxels,1)];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                      sparse(cst{i,4},1:numel(cst{i,4}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}))];
                           VoxelID                 = [VoxelID ;cst{i,4}];
                           ConstraintID            = [ConstraintID, repmat(1 + ConstraintID(end),1,numel(cst{i,4}))];

                        elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                           delta = jacobVec./(2*dij.bx(cst{i,4}).*ScaledEffect(cst{i,4}));

                           mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4},1,delta,dij.numOfVoxels,1)];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                      sparse(cst{i,4},1:numel(cst{i,4}),2*delta,dij.numOfVoxels,numel(cst{i,4}))];
                           VoxelID                 = [VoxelID ;cst{i,4}];
                           ConstraintID            = [ConstraintID, repmat(1 + ConstraintID(end),1,numel(cst{i,4}))];

                        end
                                              
                    end

                end

            end

        end

    end

end

% Calculate jacobian with dij projections
for i = 1:dij.numOfScenarios
    if isequal(type,'none')

        if ~isempty(physicalDoseProjection)
            jacob = [jacob;physicalDoseProjection' * dij.physicalDose{i}]; 
        end

    elseif isequal(type,'effect') || isequal(type,'RBExD')

        if ~isempty(mSqrtBetaDoseProjection) && ~isempty(mAlphaDoseProjection)
            mSqrtBetaDoseProjection = mSqrtBetaDoseProjection' * dij.mSqrtBetaDose{i} * w;
            mSqrtBetaDoseProjection = sparse(voxelID,constraintID(2:end),mSqrtBetaDoseProjection,...
                                         size(mAlphaDoseProjection,1),size(mAlphaDoseProjection,2));
            jacob                   = [jacob;mAlphaDoseProjection' * dij.mAlphaDose{i} + mSqrtBetaDoseProjection' * dij.mSqrtBetaDose{i}];
        end
    end
end

end