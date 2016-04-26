function g = matRad_gradFunc(w,dij,cst,type)
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
numVoxels = size(dij.physicalDose,1);

% Initializes delta
delta_underdose = zeros(numVoxels,1);
delta_overdose  = zeros(numVoxels,1);
delta_deviation = zeros(numVoxels,1);
delta_mean      = zeros(numVoxels,1);
delta_EUD       = zeros(numVoxels,1);
delta_DVH       = zeros(numVoxels,1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
 
        % get dose vector in current VOI
        d_i = d(cst{i,4});
                
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
             end
             
            if isequal(cst{i,6}(j).type, 'square underdosing') 
               
                if ~isequal(cst{i,3},'OAR')
                    % underdose : Dose minus prefered dose
                    underdose = d_i - d_ref;

                    % apply positive operator
                    underdose(underdose>0) = 0;

                    % calculate delta
                    delta_underdose(cst{i,4}) = delta_underdose(cst{i,4}) +...
                        (cst{i,6}(j).penalty/size(cst{i,4},1))*underdose;
                else
                    disp(['square underdosing constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                    % overdose : Dose minus prefered dose
                    overdose = d_i - d_ref;

                    % apply positive operator
                    overdose(overdose<0) = 0;

                    %calculate delta
                    delta_overdose(cst{i,4}) = delta_overdose(cst{i,4}) + ...
                        (cst{i,6}(j).penalty/size(cst{i,4},1))*overdose;
                
           elseif isequal(cst{i,6}(j).type, 'square deviation')
               
               if ~isequal(cst{i,3},'OAR')
                    % deviation : Dose minus prefered dose
                    deviation = d_i - d_ref;

                    % calculate delta
                    delta_deviation(cst{i,4}) = delta_deviation(cst{i,4}) +...
                        (cst{i,6}(j).penalty/size(cst{i,4},1))*deviation;
                else
                    disp(['square deviation constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                if ~isequal(cst{i,3},'TARGET')

                    % calculate delta
                    delta_mean(cst{i,4}) = delta_mean(cst{i,4}) + ...
                        (cst{i,6}(j).penalty/size(cst{i,4},1))*ones(size(cst{i,4},1),1);
                else
                    disp(['mean constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'EUD') 
               
               if ~isequal(cst{i,3},'TARGET')
                    % get exponent for EUD
                    exponent = cst{i,6}(j).EUD;

                    % calculate objective function and delta
                    if sum(d_i.^exponent)>0

                        delta_EUD(cst{i,4}) = delta_EUD(cst{i,4}) + ...
                            cst{i,6}(j).penalty*nthroot(1/size(cst{i,4},1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * (d_i.^(exponent-1));

                    end

                    if sum(~isfinite(delta_EUD)) > 0 % check for inf and nan for numerical stability
                        error(['EUD computation for ' cst{i,2} ' failed. Reduce exponent to resolve numerical problems.']);
                    end
               else
                    disp(['EUD constraint for ' cst{i,2} ' will be skipped'])
               end
                
            elseif isequal(cst{i,6}(j).type, 'max DVH objective') ||...
                   isequal(cst{i,6}(j).type, 'min DVH objective')
                
                % get reference Volume
                refVol = cst{i,6}(j).volume/100;
                
                % calc deviation
                deviation = d_i - d_ref;
                
                % calc d_ref2: V(d_ref2) = refVol
                d_ref2 = matRad_calcInversDVH(refVol,d_i);
                
                % apply lower and upper dose limits
                if isequal(cst{i,6}(j).type, 'max DVH objective')
                     deviation(d_i < d_ref | d_i > d_ref2) = 0;
                elseif isequal(cst{i,6}(j).type, 'min DVH objective')
                     deviation(d_i > d_ref | d_i < d_ref2) = 0;
                end

                % calculate delta
                delta_DVH(cst{i,4}) = delta_DVH(cst{i,4}) + (cst{i,6}(j).penalty/size(cst{i,4},1))*deviation;
                
            end
      
        end
        
    end
end
   
% Calculate gradient
if isequal(type,'none')
    
    g = (( 2*(delta_underdose + delta_overdose + delta_deviation + delta_DVH) + delta_mean + delta_EUD)' * dij.physicalDose)';

elseif isequal(type,'effect')
    
    delta = 2*(delta_underdose + delta_overdose + delta_deviation + delta_DVH) + delta_mean + delta_EUD;        
    vBias = (delta' * dij.mAlphaDose)';
    quadTerm = dij.mSqrtBetaDose * w;
    mPsi = (2*(delta.*quadTerm)'*dij.mSqrtBetaDose)';
    g    =  vBias+mPsi ; 
    
elseif isequal(type,'RBExD')
    
    ScaledEffect = d + dij.gamma;
        
    delta = 2*(delta_underdose + delta_overdose + delta_deviation + delta_DVH) + delta_mean + delta_EUD;        
    delta = delta./(2*dij.bx.*ScaledEffect);
    vBias = (delta' * dij.mAlphaDose)';
    quadTerm = dij.mSqrtBetaDose * w;
    mPsi = (2*(delta.*quadTerm)'*dij.mSqrtBetaDose)';
    g    =  vBias+mPsi ;
    
end

end