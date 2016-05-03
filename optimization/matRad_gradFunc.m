function delta = matRad_gradFunc(d_i,objective,d_ref)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: gradient function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage,
% squared deviation and DVH objectives
% 
% call
%   delta = matRad_gradFunc(d_i,objective,d_ref)
%
% input
%   d_i:       dose vector
%   objective: matRad objective struct
%   d_ref:     reference dose
%
% output
%   delta: gradient of objective function with respect to dose! needs
%   subsequent differentation for gradient in beamlet weights (see
%   gradFuncWrapper!)
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
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

numOfVoxels = numel(d_i);

if isequal(objective.type, 'square underdosing') 

    % underdose : Dose minus prefered dose
    underdose = d_i - d_ref;

    % apply positive operator
    underdose(underdose>0) = 0;

    % calculate delta
    delta = 2 * (objective.penalty/numOfVoxels)*underdose;

elseif isequal(objective.type, 'square overdosing')

    % overdose : Dose minus prefered dose
    overdose = d_i - d_ref;

    % apply positive operator
    overdose(overdose<0) = 0;

    %calculate delta
    delta = 2 * (objective.penalty/numOfVoxels)*overdose;

elseif isequal(objective.type, 'square deviation')

    % deviation : Dose minus prefered dose
    deviation = d_i - d_ref;

    % calculate delta
    delta = 2 * (objective.penalty/numOfVoxels)*deviation;

elseif isequal(objective.type, 'mean')              

    % calculate delta
    delta = (objective.penalty/numOfVoxels)*ones(numOfVoxels,1);

elseif isequal(objective.type, 'EUD') 

    % get exponent for EUD
    exponent = objective.EUD;

    % calculate objective function and delta
    if sum(d_i.^exponent)>0

        delta = objective.penalty*nthroot(1/numOfVoxels,exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * (d_i.^(exponent-1));
        
    else
        
        delta = zeros(size(d_i));

    end

    if any(~isfinite(delta)) % check for inf and nan for numerical stability
        error(['EUD computation failed. Reduce exponent to resolve numerical problems.']);
    end

elseif isequal(objective.type, 'max DVH objective') ||...
       isequal(objective.type, 'min DVH objective')

    % get reference Volume
    refVol = objective.volume/100;

    % calc deviation
    deviation = d_i - d_ref;

    % calc d_ref2: V(d_ref2) = refVol
    d_ref2 = matRad_calcInversDVH(refVol,d_i);

    % apply lower and upper dose limits
    if isequal(objective.type, 'max DVH objective')
         deviation(d_i < d_ref | d_i > d_ref2) = 0;
    elseif isequal(objective.type, 'min DVH objective')
         deviation(d_i > d_ref | d_i < d_ref2) = 0;
    end

    % calculate delta
    delta = 2 * (objective.penalty/numOfVoxels)*deviation;
    
end
     