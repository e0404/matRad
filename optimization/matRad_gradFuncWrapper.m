function g = matRad_gradFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: gradient function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage,
% squared deviation and DVH objectives
% 
% call
%   g = matRad_gradFunc(w,dij,cst,type)
%
% input
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%   type: type of optimizaiton; either 'none','effect' or 'RBExD'
%
% output
%   g: gradient of objective function
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

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Numbers of voxels
numVoxels = dij.numOfVoxels;

% Initializes delta
delta = zeros(numVoxels,1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
        
            % reference dose/effect/RBExDose
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(type,'effect') 

                    d_ref = dij.ax(cst{i,4}).*cst{i,6}(j).dose + dij.bx(cst{i,4})*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
            else
                d_ref = [];
            end
       
            if strcmp(cst{i,6}(j).robustness,'none')
                
                d_i = d{1}(cst{i,4});
                
                delta(cst{i,4}) = delta(cst{i,4}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);
                
            elseif strcmp(cst{i,6}(j).robustness,'probabilistic')
            elseif strcmp(cst{i,6}(j).robustness,'objective-wise worst case')
            elseif strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')
            end
        end
            
    end
    
end

% delta =  2*(delta_underdose + delta_overdose + delta_deviation + delta_DVH) + delta_mean + delta_EUD
  
% Calculate gradient
if isequal(type,'none')
    
    g = delta' * dij.physicalDose{1};

elseif isequal(type,'effect')
    
    vBias = (delta' * dij.mAlphaDose{1})';
    quadTerm = dij.mSqrtBetaDose{1} * w;
    mPsi = (2*(delta.*quadTerm)'*dij.mSqrtBetaDose{1})';
    g    =  vBias+mPsi ; 
    
elseif isequal(type,'RBExD')
    
    ScaledEffect = d{1} + dij.gamma;
        
    delta = delta./(2*dij.bx.*ScaledEffect);
    vBias = (delta' * dij.mAlphaDose{1})';
    quadTerm = dij.mSqrtBetaDose{1} * w;
    mPsi = (2*(delta.*quadTerm)'*dij.mSqrtBetaDose{1})';
    g    =  vBias+mPsi ;
    
end

end