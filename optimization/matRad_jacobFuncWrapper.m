function jacob = matRad_jacobFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, max DVH constraint, 
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

global matRad_DCH_ScenarioFlag;

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
scenID                  = [];
scenID2                 = [];
covConstraintID         = 0;

% compute objective function for every VOI.
for i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})

            % only perform computations for constraints
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))
                
                % compute reference
                if (~isequal(cst{i,6}(j).type, 'max dose constraint') && ~isequal(cst{i,6}(j).type, 'min dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min mean dose constraint') && ~isequal(cst{i,6}(j).type, 'max mean dose constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'min max mean dose constraint') && ~isequal(cst{i,6}(j).type, 'min EUD constraint') &&...
                    ~isequal(cst{i,6}(j).type, 'max EUD constraint') && ~isequal(cst{i,6}(j).type, 'min max EUD constraint')) &&...
                    isequal(type,'effect')
                     
                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % if conventional opt: just add constraints of nominal dose
                if strcmp(cst{i,6}(j).robustness,'none')

                    d_i = d{1}(cst{i,4}{1});

                    jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                    
                    scenID          = [scenID;1];
                    scenID2         = [scenID2;ones(numel(cst{i,4}{1}),1)];
                    covConstraintID = [covConstraintID;covConstraintID(end) + 1];
                    
                    if isequal(type,'none') && ~isempty(jacobVec)

                       physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                    elseif isequal(type,'effect') && ~isempty(jacobVec)

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                       voxelID                 = [voxelID ;cst{i,4}{1}];
                       constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    elseif isequal(type,'RBExD') && ~isempty(jacobVec)
                                        
                       scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i).^2;

                       delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);

                       mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                       mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                  sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                       voxelID                 = [voxelID ;cst{i,4}{1}];
                       constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                    end

                % if prob opt or voxel-wise worst case: add constraints of all dose scenarios
                elseif strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
                    
                    for k = 1:dij.numOfScenarios
                        
                        d_i = d{k}(cst{i,4}{1});
                        
                        jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                        
                        scenID          = [scenID;k];
                        scenID2         = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];
                        covConstraintID = [covConstraintID;covConstraintID(end) + 1];
                        
                        if isequal(type,'none') && ~isempty(jacobVec)

                           physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

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
                    
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
  
                    if isequal(cst{i,6}(j).type, 'max DCH constraint') || ...
                       isequal(cst{i,6}(j).type, 'min DCH constraint')
                    
                        % calculate scaling
                        for k = 1:dij.numOfScenarios

                            % get current dose
                            d_i = d{k}(cst{i,4}{1});

                            % inverse DVH calculation
                            d_pi(k) = matRad_calcInversDVH(cst{i,6}(j).volume/100,d_i);
                        end

                        % sort doses
                        d_pi_sort = sort(d_pi);

                        % calculate scaling
                        voxelRatio   = 1;
                        NoVoxels     = max(voxelRatio*numel(d_pi),10);
                        absDiffsort  = sort(abs(d_ref - d_pi_sort));
                        deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

                        % calclulate DVHC scaling
                        referenceVal = 0.01;
                        scaling      = min((log(1/referenceVal-1))/(2*deltaDoseMax),250);  

                        covConstraintID = [covConstraintID;repmat(1 + covConstraintID(end),dij.numOfScenarios,1)];
                                      
                        for k = 1:dij.numOfScenarios

                            d_i = d{k}(cst{i,4}{1});

                            jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,d_pi(k),scaling)/dij.numOfScenarios;

                            scenID  = [scenID;k];
                            scenID2 = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];

                            if isequal(type,'none') && ~isempty(jacobVec)

                               physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

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
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint2') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint2')         
                        
                        d_i = [];
                        
                        % get cst index of VOI that corresponds to VOI ring
                        cstidx = find(strcmp(cst(:,2),cst{i,2}(1:end-4)));
                       
                        % get dose of VOI that corresponds to VOI ring
                        for k = 1:dij.numOfScenarios
                            d_i{k} = d{k}(cst{cstidx,4}{1});
                        end

                        % calc invers DCH of VOI
                        refQ   = cst{i,6}(j).coverage/100;
                        refVol = cst{i,6}(j).volume/100;
                        d_ref2 = matRad_calcInversDCH(refVol,refQ,d_i,dij.numOfScenarios);

                        % get dose of Target Ring
                        d_i = d{1}(cst{i,4}{1});

                        % calc voxel dependent weighting
                        voxelWeighting = 5*cst{i,5}.voxelProb;

                        jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref,1,1,d_ref2,voxelWeighting);

                        scenID          = [scenID;1];
                        scenID2         = [scenID2;ones(numel(cst{i,4}{1}),1)];
                        covConstraintID = [covConstraintID;covConstraintID(end) + 1];

                        if isequal(type,'none') && ~isempty(jacobVec)

                           physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

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
                        
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint3') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint3')
                                              
                        % calculate scenario approximation scaling
                        for k = 1:dij.numOfScenarios

                            % get current dose
                            d_i = d{k}(cst{i,4}{1});

                            % calculate volume
                            volume_pi(k) = sum(d_i >= d_ref)/numel(d_i);
                        end
                        
                        % sort columes
                        volumes_pi_sort = sort(volume_pi);

                        % calculate scaling
                        voxelRatio   = cst{i,6}(j).coverage/100;
                        NoVoxels     = voxelRatio*numel(volume_pi);
                        absDiffsort  = sort(abs(cst{i,6}(j).volume/100 - volumes_pi_sort));
                        deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

                        % calclulate DVHC scaling
                        referenceVal = 0.01;
                        scaling      = min((log(1/referenceVal-1))/(2*deltaDoseMax),250);
                        
                        scenProb        = 1/dij.numOfScenarios;     % assume scenarios with equal probabilities
                        covConstraintID = [covConstraintID;repmat(1 + covConstraintID(end),dij.numOfScenarios,1)];
                       
                        for k = 1:dij.numOfScenarios
                            
                            d_i = d{k}(cst{i,4}{1});
                            
                            jacobVec = scenProb*2*scaling*exp(2*scaling*(volume_pi(k)-cst{i,6}(j).volume/100))./(exp(2*scaling*(volume_pi(k)-cst{i,6}(j).volume/100))+1).^2;

                            jacobVec =  jacobVec*matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);

%                             jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
%                             covConstraintID = [covConstraintID;covConstraintID(end) + 1];

                            scenID  = [scenID;k];
                            scenID2 = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];

                            if isequal(type,'none') && ~isempty(jacobVec)

                               physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                            elseif isequal(type,'effect') && ~isempty(jacobVec)

                               mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                               mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                          sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                               voxelID                 = [voxelID ;cst{i,4}{1}];
                               constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                            elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                               scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i).^2;

                               delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);

                               mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                               mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                          sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                               voxelID                 = [voxelID ;cst{i,4}{1}];
                               constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                            end
                            
                        end
                        
                    elseif isequal(cst{i,6}(j).type, 'max DCH constraint4') || ...
                           isequal(cst{i,6}(j).type, 'min DCH constraint4')
                                              
                        for k = 1:dij.numOfScenarios
                            
                            d_i = d{k}(cst{i,4}{1});
                            
                            if matRad_DCH_ScenarioFlag(k)
                                jacobVec =  matRad_jacobFunc(d_i,cst{i,6}(j),d_ref);
                            else
                                jacobVec = 0;
                            end
                            
                            covConstraintID = [covConstraintID;covConstraintID(end) + 1];

                            scenID  = [scenID;k];
                            scenID2 = [scenID2;repmat(k,numel(cst{i,4}{1}),1)];

                            if isequal(type,'none') && ~isempty(jacobVec)

                               physicalDoseProjection = [physicalDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];

                            elseif isequal(type,'effect') && ~isempty(jacobVec)

                               mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,jacobVec,dij.numOfVoxels,1)];
                               mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                          sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*jacobVec,dij.numOfVoxels,numel(cst{i,4}{1}))];
                               voxelID                 = [voxelID ;cst{i,4}{1}];
                               constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                            elseif isequal(type,'RBExD') && ~isempty(jacobVec)

                               scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i).^2;

                               delta = jacobVec./(2*dij.bx(cst{i,4}{1}).*scaledEffect);

                               mAlphaDoseProjection    = [mAlphaDoseProjection,sparse(cst{i,4}{1},1,delta,dij.numOfVoxels,1)];
                               mSqrtBetaDoseProjection = [mSqrtBetaDoseProjection,...
                                                          sparse(cst{i,4}{1},1:numel(cst{i,4}{1}),2*delta,dij.numOfVoxels,numel(cst{i,4}{1}))];
                               voxelID                 = [voxelID ;cst{i,4}{1}];
                               constraintID            = [constraintID, repmat(1 + constraintID(end),1,numel(cst{i,4}{1}))];

                            end
                            
                        end                        
                       
                    end

                end

            end

        end

    end

end

if isequal(type,'effect') || isequal(type,'RBExD')
    constraintID = constraintID(2:end);
end

% Calculate jacobian with dij projections
for i = 1:dij.numOfScenarios
    if isequal(type,'none')

        if ~isempty(physicalDoseProjection)
            
            jacobLogical          = (scenID == i);
            jacob(jacobLogical,:) = physicalDoseProjection(:,jacobLogical)' * dij.physicalDose{i};
            
        end

    elseif isequal(type,'effect') || isequal(type,'RBExD')

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

% summ over coverage scenarios
if numel(unique(covConstraintID)) < numel(covConstraintID)
    
    [rows, cols] = ndgrid(covConstraintID(2:end),1:size(jacob,2));
    jacob        = accumarray([rows(:) cols(:)],jacob(:));
   
end



end
