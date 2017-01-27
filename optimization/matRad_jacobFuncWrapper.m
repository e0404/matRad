function jacob = matRad_jacobFuncWrapper(w,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min mean dose constraint, max mean dose constraint,
% min EUD constraint, max EUD constraint, max DVH constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(w,dij,cst,options)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   options: option struct defining the type of optimization
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

global cScaling
global kDVH;
global kDCH;
global matRad_iteration;
global JACOBIAN;

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,options);

% initialize jacobian
jacob = sparse([]);

% initialize projection matrices and id containers
DoseProjection          = sparse([]);
mAlphaDoseProjection    = sparse([]);
mSqrtBetaDoseProjection = sparse([]);
voxelID                 = [];
constraintID            = 0;
constraintID2           = 0;
scenID                  = [];
scenID2                 = [];

% compute objective function for every VOI.
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'max dose constraint')      && ~isequal(cst{i,6}(j).type, 'min dose constraint')          &&...
                    ~isequal(cst{i,6}(j).type, 'max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ...
                    ~isequal(cst{i,6}(j).type, 'min EUD constraint')       && ~isequal(cst{i,6}(j).type, 'max EUD constraint'))           && ...
                    ( isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LSM_effect'))
                     
                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                    
                    % set ID vectors
                    scenID        = [scenID;1];
                    scenID2       = [scenID2;ones(numel(cst{i,4}{1}),1)];
                    constraintID2 = [constraintID2;constraintID2(end) + 1];
                    
                    if isequal(options.bioOpt,'none') && ~isempty(jacobVec) || isequal(options.ID,'protons_const_RBExD')

                       DoseProjection          = [DoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                    elseif (isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LSM_effect')) && ~isempty(jacobVec)

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                       voxelID                 = [voxelID ;cst{i,4}{1}];
                       constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    elseif (isequal(options.bioOpt,'LEMIV_RBExD') || isequal(options.bioOpt,'LSM_RBExD')) && ~isempty(jacobVec)
                                        
                       scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);

                       delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                       voxelID                 = [voxelID ;cst{i,4}{1}];
                       constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    end

                % if prob opt or voxel-wise worst case: add constraints of all dose scenarios
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'WC')
                    
                    for k = 1:dij.numOfScenarios
                        
                        d_i = d{k}(cst{i,4}{1});
                        
                        jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                        
                        % set ID vectors
                        scenID        = [scenID;k];
                        scenID2       = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];
                        constraintID2 = [constraintID2;constraintID2(end) + 1];
                        
                        if isequal(type,'none') && ~isempty(jacobVec)

                           DoseProjection = [DoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                        elseif isequal(type,'effect') && ~isempty(jacobVec)

                           mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                      sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                           voxelID                 = [voxelID ;cst{i,4}{1}];
                           constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                        elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                           delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*ScaledEffect(cst{i,4}{1}));

                           mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                           mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                      sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                           voxelID                 = [voxelID ;cst{i,4}{1}];
                           constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                        end
                                              
                    end
                    
                % if coveraged based opt     
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
  
                    if isequal(cst{i,6}(j).type, 'max DCH Area constraint') || ...
                       isequal(cst{i,6}(j).type, 'min DCH Area constraint')
                   
                        % calc invers DCH
                        Q_ref  = cst{i,6}(j).coverage/100;
                        V_ref  = cst{i,6}(j).volume/100;
                        d_ref2 = matRad_calcInversDCH(V_ref,Q_ref,d,dij,cst(i,:));
                            
                        if dij.numOfScenarios > 1
                            
                            % set coverage constraint ID2 vector
                            constraintID2 = [constraintID2;repmat(1 + constraintID2(end),dij.numOfScenarios,1)];

                            for k = 1:dij.numOfScenarios

                                % get VOI dose in current scenario
                                d_i = d{k}(cst{i,4}{1});

                                % get voxel dependent weigthing
                                voxelWeighting = 1;  

                                % calculate delta
                                jacobVec = dij.ScenProb(k) * matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);                   
                                
                                % set ID vectors
                                scenID  = [scenID;k];
                                scenID2 = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];
                                
                                if isequal(type,'none') && ~isempty(jacobVec)

                                   DoseProjection = [DoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                                elseif isequal(type,'effect') && ~isempty(jacobVec)

                                    error('effect optimization in Area constraint not implemented yet')

                                elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                                    error('RBExD optimization in Area constraint not implemented yet')

                                end                                

                            end

                        else                                                    

                            % get VOI ScenUnion dose of nominal scneario
                            cstLogical = strcmp(cst(:,2),[cst{i,2},' ScenUnion']);
                            d_i        = d{1}(cst{cstLogical,5}.voxelID);

                            % get voxel dependent weigthing
                            voxelWeighting = 1;  

                            % calculate delta
                            jacobVec = matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);                  
                            
                            % set ID vectors
                            scenID        = [scenID;1];
                            scenID2       = [scenID2;ones(numel(cst{i,4}{1}),1)];
                            constraintID2 = [constraintID2;constraintID2(end) + 1];

                            if isequal(type,'none') && ~isempty(jacobVec)

                               DoseProjection = [DoseProjection,sparse(cst{cstLogical,5}.voxelID,1,jacobVec,dij.numOfVoxels,1)];

                            elseif isequal(type,'effect') && ~isempty(jacobVec)

                                error('effect optimization in Area constraint not implemented yet')

                            elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                                error('RBExD optimization in Area constraint not implemented yet')

                            end

                        end
                        
                    elseif isequal(cst{i,6}(j).type, 'max DCH Theta constraint') || ...
                           isequal(cst{i,6}(j).type, 'min DCH Theta constraint')
                                              
                        if dij.numOfScenarios > 1
                            
                            % get VOI dose union in every scneario
                            dUnion = [];
                            for k = 1:dij.numOfScenarios
                                % get current dose
                                dUnion = union(dUnion,d{k}(cst{i,4}{1}));
                            end
                            
                            % calculate DVH logistic function scaling
                            DVHScaling = matRad_calcLogisticFuncScaling(dUnion,d_ref,0.5,0.01,0,250);                            
                            kDVH(j,matRad_iteration+1)= DVHScaling;
                            
                            for k = 1:dij.numOfScenarios

                                % get VOI dose in current scenario
                                d_i = d{k}(cst{i,4}{1});

                                % calculate volumes
                                volume_pi(k) = sum(1./(1+exp(-2*DVHScaling*(d_i-d_ref))))/numel(d_i);

                            end
                            
                            % get scenario probabilities
                            scenProb = dij.ScenProb;
                            
                        else
                            
                            % get VOI ScenUnion dose of nominal scneario
                            cstLogical = strcmp(cst(:,2),[cst{i,2},' ScenUnion']);
                            d_i        = d{1}(cst{cstLogical,5}.voxelID);
                            
                            % calculate DVH logistic function scaling
                            DVHScaling = matRad_calcLogisticFuncScaling(d_i,d_ref,0.5,0.01,0,250);                            
                            kDVH(j,matRad_iteration+1)= DVHScaling;
                            
                            for k = 1:cst{i,5}.VOIShift.ncase

                                % get VOI dose in current scenario
                                if isequal(cst{i,5}.VOIShift.shiftType,'rounded')
                                    d_i = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));
                                elseif isequal(cst{i,5}.VOIShift.shiftType,'linInterp')
                                    error('linInterp in DCH Theta constraint not implemented yet')
                                end
                                
                                % calculate volumes
                                volume_pi(k) = sum(1./(1+exp(-2*DVHScaling*(d_i-d_ref))))/numel(d_i);

                            end
                            
                            % get scenario probabilities
                            scenProb = 1/cst{i,5}.VOIShift.ncase;  % assume equiprobable scenarios
                            
                        end

                        % calculate DCH logistic function scaling
                        DCHScaling = matRad_calcLogisticFuncScaling(volume_pi,cst{i,6}(j).volume/100,1,0.01,0,1000);
                        kDCH(j,matRad_iteration+1)= DCHScaling;
                        
                        if dij.numOfScenarios > 1
                            
                            % set coverage constraint ID2 vector
                            constraintID2 = [constraintID2;repmat(1 + constraintID2(end),dij.numOfScenarios,1)];
                            
                            for k = 1:dij.numOfScenarios
                                
                                % get VOI dose in current scenario
                                d_i = d{k}(cst{i,4}{1});
                                
                                % calculate delta
                                jacobVec = scenProb(k)*matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,0,0,DVHScaling,DCHScaling,volume_pi(k));
                                    
                                % set ID vectors
                                scenID  = [scenID;k];
                                scenID2 = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];

                                if isequal(type,'none') && ~isempty(jacobVec)

                                   DoseProjection = [DoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                                elseif isequal(type,'effect') && ~isempty(jacobVec)

                                    error('effect optimization in Theta constraint not implemented yet')

                                elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                                    error('RBExD optimization in Theta constraint not implemented yet')

                                end

                            end
                        else
                            
                            % set coverage constraint ID2 vector
                            constraintID2 = [constraintID2;repmat(1 + constraintID2(end),cst{i,5}.VOIShift.ncase,1)];
                            
                            for k = 1:cst{i,5}.VOIShift.ncase
                                
                                % get VOI dose in current scenario
                                d_i = d{1}(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k));
                                
                                % calculate delta
                                jacobVec = scenProb*matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,0,0,DVHScaling,DCHScaling,volume_pi(k));
                                    
                                % set ID vectors
                                scenID  = [scenID;1];
                                scenID2 = [scenID2;ones(numel(cst{i,4}{1}),1)];

                                if isequal(type,'none') && ~isempty(jacobVec)

                                   DoseProjection = [DoseProjection,sparse(cst{i,4}{1}-cst{i,5}.VOIShift.roundedShift.idxShift(k),1,jacobVec,dij.numOfVoxels,1)];

                                elseif isequal(type,'effect') && ~isempty(jacobVec)

                                    error('effect optimization in Theta constraint not implemented yet')

                                elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                                    error('RBExD optimization in Theta constraint not implemented yet')

                                end

                            end
                        end
                        
                    end

                end

            end

        end

    end

end

if isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD')
    constraintID = constraintID(2:end);
end

% Calculate jacobian with dij projections
for i = 1:dij.numOfScenarios
   % enter if statement also for protons using a constant RBE
   if isequal(options.bioOpt,'none') ||  isequal(options.ID,'protons_const_RBExD')

        if ~isempty(DoseProjection)
            
            jacobLogical          = (scenID == i);
            jacob(jacobLogical,:) = DoseProjection(:,jacobLogical)' * dij.physicalDose{i};
            
        end

    elseif isequal(options.bioOpt,'LEMIV_effect') || isequal(options.bioOpt,'LEMIV_RBExD')

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

% summ over scenarios with eqaul constraintID2
if numel(unique(constraintID2)) < numel(constraintID2)
    
    [rows, cols] = ndgrid(constraintID2(2:end),1:size(jacob,2));
    jacob        = accumarray([rows(:) cols(:)],jacob(:));
   
end

% apply constraint scaling
jacob = bsxfun(@times,cScaling,jacob);

% save min/max jacobian
if ~isempty(jacob)
    for k = 1:size(jacob,1)
        jacob_ = jacob(k,:);
        if isempty(max(abs(jacob_(jacob_ ~= 0))))
            JACOBIAN(k,1,matRad_iteration+1)= 0;
        else
            JACOBIAN(k,1,matRad_iteration+1)= max(abs(jacob_(jacob_ ~= 0)));
        end
        if isempty(min(abs(jacob_(jacob_ ~= 0))))
            JACOBIAN(k,2,matRad_iteration+1)= 0;
        else
            JACOBIAN(k,2,matRad_iteration+1)= min(abs(jacob_(jacob_ ~= 0)));
        end
    end
end

end
