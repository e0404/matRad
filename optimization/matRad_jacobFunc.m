function jacobVec = matRad_jacobFunc(d_i,constraint,d_ref,d_pi,scaling,d_ref2,weighting)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: jacobian function for inverse planning supporting max dose
% constraint, min dose constraint, min max dose constraint, min mean, max
% min, min max mean constraint, min EUD constraint, max EUDconstraint, 
% min max EUD constraint, max DVH constraint, 
% min DVH constraint 
% 
% call
%   jacobVec = matRad_jacobFunc(d_i,constraint,d_ref)
%
% input
%   d_i:        dose vector
%   constraint: matRad constraint struct
%   d_ref:      reference dose
%
% output
%   jacobVec:  jacobian vector of constraint for differentation with
%              respect to dose. need subsequent differentation for jacobian
%              in beamlet weights (see jacobFunWrapper)
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

global matRad_DVH_Scaling;

numOfVoxels = numel(d_i);

if isequal(constraint.type, 'max dose constraint')
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837

    epsilon = 1e-3;

    jacobVec = exp( (d_i-max(d_i))/epsilon );
    jacobVec = jacobVec/sum(jacobVec);

elseif isequal(constraint.type, 'min dose constraint')
    % use log sum exp approximation, see appendix A in
    % http://scitation.aip.org/content/aapm/journal/medphys/41/8/10.1118/1.4883837

    epsilon = 1e-3;

    jacobVec = exp( (min(d_i)-d_i)/epsilon );
    jacobVec = jacobVec/sum(jacobVec);

elseif isequal(constraint.type, 'max mean dose constraint') || ...
       isequal(constraint.type, 'min mean dose constraint') || ...
       isequal(constraint.type, 'min max mean dose constraint') 

    jacobVec = ones(numOfVoxels,1)./numOfVoxels;

elseif isequal(constraint.type, 'max EUD constraint') || ...
       isequal(constraint.type, 'min EUD constraint') || ...
       isequal(constraint.type, 'min max EUD constraint') 

    % exponenent for EUD constraint
    exponent = constraint.EUD;

    jacobVec = nthroot(1/numOfVoxels,exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * ...
                  (d_i.^(exponent-1));

elseif isequal(constraint.type, 'max DVH constraint') || ...
       isequal(constraint.type, 'min DVH constraint')

    d_i_sort = sort(d_i);

    % calculate scaling
    VoxelRatio   = 1;
    NoVoxels     = max(VoxelRatio*numel(d_i),10);
    absDiffsort  = sort(abs(d_ref - d_i_sort));
    deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

    % calclulate DVHC scaling
    ReferenceVal            = 0.01;
    DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);

    jacobVec = (2/numOfVoxels)*DVHCScaling*exp(2*DVHCScaling*(d_i-d_ref))./(exp(2*DVHCScaling*(d_i-d_ref))+1).^2;

    % alternative constraint calculation 4/4 %               
    % % get reference Volume
    % refVol = cst{j,6}(k).volume/100;
    %  
    % % calc deviation
    % deviation = d_i - d_ref;
    % 
    % % calc d_ref2: V(d_ref2) = refVol
    % d_ref2 = matRad_calcInversDVH(refVol,d_i);
    % 
    % % apply lower and upper dose limits
    % if isequal(cst{j,6}(k).type, 'max DVH constraint')
    %      deviation(d_i < d_ref | d_i > d_ref2) = 0;
    % elseif isequal(cst{j,6}(k).type, 'min DVH constraint')
    %      deviation(d_i > d_ref | d_i < d_ref2) = 0;
    % end
    %   
    % %jacobVec = ones(size(cst{j,4}));             % linear deviation
    % %jacobVec = 2*deviation;                      % square deviation
    % jacobVec = (1/size(cst{j,4},1))*2*deviation; % square deviation with normalization
    % %jacobVec = 4*(deviation).^3;                  % squared square devioation
    % alternative constraint calculation 4/4 %
    
elseif isequal(constraint.type, 'max DCH constraint') || ...
       isequal(constraint.type, 'min DCH constraint')
   
    % sort dose values
    d_i_tmp = sort(d_i,'descend');

    ix = max([1 ceil(constraint.volume/100*numel(d_i_tmp))]);
    
    d_i_tmp = d_i_tmp(ix);

    d_i(d_i > d_i_tmp) = 0;
   
    % max part
    epsilon = 1e-3;

    jacobVec = exp( (min(d_i)-d_i)/epsilon );
    jacobVec = jacobVec/sum(jacobVec);
    
    % DVH part
    jacobVec = jacobVec.*2*scaling*exp(2*scaling*(d_pi-d_ref))./(exp(2*scaling*(d_pi-d_ref))+1).^2;
    
elseif isequal(constraint.type, 'max DCH constraint2') || ...
       isequal(constraint.type, 'min DCH constraint2')
   
    % calc deviation
    deviation = d_i - d_ref;

    % apply lower and upper dose limits
    if isequal(constraint.type, 'max DCH constraint2')
     deviation(d_i < d_ref | d_i > d_ref2) = 0;
    elseif isequal(constraint.type, 'min DCH constraint2')
     deviation(d_i > d_ref | d_i < d_ref2) = 0;
    end

    % apply weighting
    deviation = deviation.*(weighting).^2';

    % calculate delta
    jacobVec = 2 * (1/numOfVoxels)*deviation;
    
elseif isequal(constraint.type, 'max DCH constraint3') || ...
       isequal(constraint.type, 'min DCH constraint3')
   
    % calculate logistic function scaling and jacobian
    DVHScaling = matRad_DVH_Scaling;    
    jacobVec   = (2/numOfVoxels)*DVHScaling*exp(2*DVHScaling*(d_i-d_ref))./(exp(2*DVHScaling*(d_i-d_ref))+1).^2; 
    
elseif isequal(constraint.type, 'max DCH constraint4') || ...
       isequal(constraint.type, 'min DCH constraint4')

    d_i_sort = sort(d_i);

    % calculate scaling
    VoxelRatio   = 1;
    NoVoxels     = max(VoxelRatio*numel(d_i),10);
    absDiffsort  = sort(abs(d_ref - d_i_sort));
    deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

    % calclulate DVHC scaling
    ReferenceVal            = 0.01;
    DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);

    jacobVec = (2/numOfVoxels)*DVHCScaling*exp(2*DVHCScaling*(d_i-d_ref))./(exp(2*DVHCScaling*(d_i-d_ref))+1).^2; 
    
elseif isequal(constraint.type, 'max DCH constraint5') || ...
       isequal(constraint.type, 'min DCH constraint5')

    d_i_sort = sort(d_i);

    % calculate scaling
    VoxelRatio   = 1;
    NoVoxels     = max(VoxelRatio*numel(d_i),10);
    absDiffsort  = sort(abs(d_ref - d_i_sort));
    deltaDoseMax = absDiffsort(ceil(NoVoxels/2));

    % calclulate DVHC scaling
    ReferenceVal            = 0.01;
    DVHCScaling             = min((log(1/ReferenceVal-1))/(2*deltaDoseMax),250);

    jacobVec = (2/numOfVoxels)*DVHCScaling*exp(2*DVHCScaling*(d_i-d_ref))./(exp(2*DVHCScaling*(d_i-d_ref))+1).^2;      

else

    jacobVec = [];

end