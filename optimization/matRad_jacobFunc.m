function jacob = matRad_jacobFunc(w,dij,cst,type)
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

% loop over all scenarios
for i = 1:dij.numOfScenarios

    % initialize projection matrices and id containers
    physicalDoseProjection  = sparse([]);
    mAlphaDoseProjection    = sparse([]);
    mSqrtBetaDoseProjection = sparse([]);
    voxelID                 = [];
    constraintID            = 0;

    % calulate Scaled effect
    if isequal(type,'RBExD')
        scaledEffect = d{i} + dij.gamma;
    end

    % compute objective function for every VOI.
    for  j = 1:size(cst,1)

        % Only take OAR or target VOI.
        if ~isempty(cst{j,4}) && ( isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') )

            % loop over the number of constraints for the current VOI
            for k = 1:numel(cst{j,6})
                
                % only in the nominal case or for robust optimization
                if i == 1 || strcmp(cst{j,6}(k).robustness,'probabilistic') || ...
                             strcmp(cst{j,6}(k).robustness,'voxel-wise worst case')

                    if isequal(cst{j,6}(k).type, 'max dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'min dose constraint') || ...
                       isequal(cst{j,6}(k).type, 'min max dose constraint') 

                       if isequal(type,'none')

                           physicalDoseProjection = [physicalDoseProjection,sparse(cst{j,4},1:numel(cst{j,4}),1,dij.numOfVoxels,numel(cst{j,4}))];

                       elseif isequal(type,'effect')

                           mAlphaDoseProjection    = [mAlphaDoseProjection,...
                                                      sparse(cst{j,4},1:numel(cst{j,4}),1,dij.numOfVoxels,numel(cst{j,4}))];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                      sparse(cst{j,4},1:numel(cst{j,4}),2,dij.numOfVoxels,numel(cst{j,4}))]; 
                           voxelID                 = [voxelID ;cst{j,4}];
                           constraintID            = [constraintID, [1:numel(cst{j,4})] + constraintID(end)];

                       elseif isequal(type,'RBExD')

                           mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1:numel(cst{j,4}),1./(2*dij.bx(cst{j,4}).*scaledEffect(cst{j,4})),...
                                                                                  dij.numOfVoxels,numel(cst{j,4}))];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,sparse(cst{j,4},1:numel(cst{j,4}),2./(2*dij.bx(cst{j,4}).*scaledEffect(cst{j,4})),...
                                                                                     dij.numOfVoxels,numel(cst{j,4}))]; 
                           voxelID                 = [voxelID ;cst{j,4}];
                           constraintID            = [constraintID, [1:numel(cst{j,4})] + constraintID(end)];

                       end

                    elseif isequal(cst{j,6}(k).type, 'max mean dose constraint') || ...
                           isequal(cst{j,6}(k).type, 'min mean dose constraint') || ...
                           isequal(cst{j,6}(k).type, 'min max mean dose constraint') 

                        if isequal(type,'none')

                            physicalDoseProjection = [physicalDoseProjection,sparse(cst{j,4},1,1/length(cst{j,4}),dij.numOfVoxels,1)];

                        elseif isequal(type,'effect')

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,1/length(cst{j,4}),dij.numOfVoxels,1)]; 
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2/numel(cst{j,4}),dij.numOfVoxels,numel(cst{j,4}))];                    
                            voxelID                 = [voxelID ;cst{j,4}];  
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        elseif isequal(type,'RBExD')

                            delta = (ones(numel(cst{j,4}),1)./numel(cst{j,4}))./(2*dij.bx(cst{j,4}).*scaledEffect(cst{j,4}));

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,delta,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2*delta,dij.numOfVoxels,numel(cst{j,4}))];
                            voxelID                 = [voxelID ;cst{j,4}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        end

                    elseif isequal(cst{j,6}(k).type, 'max EUD constraint') || ...
                           isequal(cst{j,6}(k).type, 'min EUD constraint') || ...
                           isequal(cst{j,6}(k).type, 'min max EUD constraint') 

                        % exponenent for EUD constraint
                        exponent = cst{j,6}(k).EUD;

                        % get dose in VOI
                        d_i = d{i}(cst{j,4});

                        eudJacobVec = nthroot(1/size(cst{j,4},1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * ...
                                      (d_i.^(exponent-1));

                        if isequal(type,'none')

                            physicalDoseProjection = [physicalDoseProjection,sparse(cst{j,4},1,eudJacobVec,dij.numOfVoxels,1)];

                        elseif isequal(type,'effect')

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,eudJacobVec,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2*eudJacobVec,dij.numOfVoxels,numel(cst{j,4}))];
                            voxelID                 = [voxelID ;cst{j,4}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        elseif isequal(type,'RBExD')

                            delta = eudJacobVec./(2*dij.bx(cst{j,4}).*scaledEffect(cst{j,4}));

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,delta,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2*delta,dij.numOfVoxels,numel(cst{j,4}))];
                            voxelID                 = [voxelID ;cst{j,4}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        end

                    elseif isequal(cst{j,6}(k).type, 'exact DVH constraint') || ...
                           isequal(cst{j,6}(k).type, 'max DVH constraint') || ...
                           isequal(cst{j,6}(k).type, 'min DVH constraint')

                        % get dose in VOI
                        d_i = d{i}(cst{j,4});
                        d_i_sort = sort(d_i);

                        % reference dose/effect/dosexRBE
                        if isequal(type,'effect')
                            d_ref = dij.ax(cst{j,4}).*cst{j,6}(k).dose + dij.bx(cst{j,4})*cst{j,6}(k).dose^2;
                        else
                            d_ref = cst{j,6}(k).dose;
                        end

                        % calculate scaling
                        voxelRatio   = 1;
                        noVoxels     = max(voxelRatio*numel(d_i),10);
                        absDiffsort  = sort(abs(d_ref - d_i_sort));
                        deltaDoseMax = absDiffsort(ceil(noVoxels/2));

                        % calclulate DVHC scaling
                        ReferenceVal            = 0.01;
                        DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);

                        DVHJacobVec = (2/size(cst{j,4},1))*DVHCScaling*exp(2*DVHCScaling*(d_i-d_ref))./(exp(2*DVHCScaling*(d_i-d_ref))+1).^2;

                        % alternative constraint calculation 4/4 %               
                        % % get reference Volume
                        % refVol = cst{i,6}(j).volume/100;
                        %  
                        % % calc deviation
                        % deviation = d_i - d_ref;
                        % 
                        % % calc d_ref2: V(d_ref2) = refVol
                        % d_ref2 = matRad_calcInversDVH(refVol,d_i);
                        % 
                        % % apply lower and upper dose limits
                        % if isequal(cst{i,6}(j).type, 'max DVH constraint')
                        %      deviation(d_i < d_ref | d_i > d_ref2) = 0;
                        % elseif isequal(cst{i,6}(j).type, 'min DVH constraint')
                        %      deviation(d_i > d_ref | d_i < d_ref2) = 0;
                        % end
                        %   
                        % %DVHJacobVec = ones(size(cst{i,4}));             % linear deviation
                        % %DVHJacobVec = 2*deviation;                      % square deviation
                        % DVHJacobVec = (1/size(cst{i,4},1))*2*deviation; % square deviation with normalization
                        % %DVHJacobVec = 4*(deviation).^3;                  % squared square devioation
                        % alternative constraint calculation 4/4 %

                        if isequal(type,'none')

                            physicalDoseProjection = [physicalDoseProjection,sparse(cst{j,4},1,DVHJacobVec,dij.numOfVoxels,1)];

                        elseif isequal(type,'effect')

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,DVHJacobVec,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2*DVHJacobVec,dij.numOfVoxels,numel(cst{j,4}))];
                            voxelID                 = [voxelID ;cst{j,4}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        elseif isequal(type,'RBExD')

                            delta             = DVHJacobVec./(2*dij.bx(cst{j,4}).*scaledEffect(cst{j,4}));

                            mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{j,4},1,delta,dij.numOfVoxels,1)];
                            mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                       sparse(cst{j,4},1:numel(cst{j,4}),2*delta,dij.numOfVoxels,numel(cst{j,4}))];
                            voxelID                 = [voxelID ;cst{j,4}];
                            constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{j,4}))];

                        end
                        
                    end
                    
                end

            end

        end

    end

    % Calculate jacobian with dij projections
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