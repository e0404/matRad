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

% Initializes constraints
jacob                   = sparse([]);
physicalDoseProjection  = sparse([]);
mAlphaDoseProjection    = sparse([]);
mSqrtBetaDoseProjection = sparse([]);
VoxelID                 = [];
ConstraintID            = 0;

% calulate Scaled effect
if isequal(type,'RBExD')
    ScaledEffect = d + dij.gamma;
end

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
                    
        % loop over all objectives and constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            if isequal(cst{i,6}(j).type, 'max dose constraint')
                % use log sum exp approximation, see appendix A in
                % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
                
                epsilon = 1e-3;
                   
                % get dose in VOI
                d_i = d(cst{i,4});
                
                jacobVec = exp( (d_i-max(d_i))/epsilon );
                jacobVec = jacobVec/sum(jacobVec);
                
            elseif isequal(cst{i,6}(j).type, 'min dose constraint')
                % use log sum exp approximation, see appendix A in
                % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837
                
                epsilon = 1e-3;
                   
                % get dose in VOI
                d_i = d(cst{i,4});
                
                jacobVec = exp( (min(d_i)-d_i)/epsilon );
                jacobVec = jacobVec/sum(jacobVec);
                
            elseif isequal(cst{i,6}(j).type, 'max mean dose constraint') || ...
                   isequal(cst{i,6}(j).type, 'min mean dose constraint') || ...
                   isequal(cst{i,6}(j).type, 'min max mean dose constraint') 
               
                    jacobVec = ones(numel(cst{i,4}),1)./numel(cst{i,4});
                               
            elseif isequal(cst{i,6}(j).type, 'max EUD constraint') || ...
                   isequal(cst{i,6}(j).type, 'min EUD constraint') || ...
                   isequal(cst{i,6}(j).type, 'min max EUD constraint') 
                
                % exponenent for EUD constraint
                exponent = cst{i,6}(j).EUD;

                % get dose in VOI
                d_i = d(cst{i,4});
                
                jacobVec = nthroot(1/size(cst{i,4},1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * ...
                              (d_i.^(exponent-1));
                                
            elseif isequal(cst{i,6}(j).type, 'exact DVH constraint') || ...
                   isequal(cst{i,6}(j).type, 'max DVH constraint') || ...
                   isequal(cst{i,6}(j).type, 'min DVH constraint')

                % get dose in VOI
                d_i = d(cst{i,4});
                d_i_sort = sort(d_i);
                
                % reference dose/effect/dosexRBE
                if isequal(type,'effect')
                    d_ref = dij.ax(cst{i,4}).*cst{i,6}(j).dose + dij.bx(cst{i,4})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % calculate scaling
                VoxelRatio   = 1;
                NoVoxels     = max(VoxelRatio*numel(d_i),10);
                absDiffsort  = sort(abs(d_ref - d_i_sort));
                deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

                % calclulate DVHC scaling
                ReferenceVal            = 0.01;
                DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);
                
                jacobVec = (2/size(cst{i,4},1))*DVHCScaling*exp(2*DVHCScaling*(d_i-d_ref))./(exp(2*DVHCScaling*(d_i-d_ref))+1).^2;
                
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
                % %jacobVec = ones(size(cst{i,4}));             % linear deviation
                % %jacobVec = 2*deviation;                      % square deviation
                % jacobVec = (1/size(cst{i,4},1))*2*deviation; % square deviation with normalization
                % %jacobVec = 4*(deviation).^3;                  % squared square devioation
                % alternative constraint calculation 4/4 %
            
            else
            
                jacobVec = [];
                
            end
                                          
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

% Calculate jacobian with dij projections
if isequal(type,'none')
    
    if isempty(physicalDoseProjection)
        jacob = [];
    else   
        jacob = physicalDoseProjection' * dij.physicalDose; 
    end

elseif isequal(type,'effect') || isequal(type,'RBExD')
    
    if isempty(mSqrtBetaDoseProjection) || isempty(mAlphaDoseProjection)
        jacob = [];
    else
    mSqrtBetaDoseProjection = mSqrtBetaDoseProjection' * dij.mSqrtBetaDose * w;
    mSqrtBetaDoseProjection = sparse(VoxelID,ConstraintID(2:end),mSqrtBetaDoseProjection,...
                                     size(mAlphaDoseProjection,1),size(mAlphaDoseProjection,2));
    jacob                   = mAlphaDoseProjection' * dij.mAlphaDose + mSqrtBetaDoseProjection' * dij.mSqrtBetaDose;
    end
end

end