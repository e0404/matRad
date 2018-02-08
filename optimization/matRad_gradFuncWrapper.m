function g = matRad_gradFuncWrapper(w,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: gradient function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage,
% squared deviation and DVH objectives
% 
% call
%   g = matRad_gradFuncWrapper(w,dij,cst,options)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   options: option struct defining the type of optimization
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

% read in global dose and bixel variables
global matRad_global_x;
global matRad_global_d;

if ~isequal(w,matRad_global_x)
    % new bixel weights, update dose
    % get current dose / effect / RBExDose vector
    d = matRad_backProjection(w,dij,options);
    matRad_global_d = d;
    matRad_global_x = w;
else
    % old bixel weights, use global dose
    d = matRad_global_d;
end

% Initializes delta
delta      = cell(options.numOfScenarios,1);
[delta{:}] = deal(zeros(dij.numOfVoxels,1));

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
        
            % only perform gradient computations for objectives
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                % compute reference
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(options.bioOpt,'LEMIV_effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % different gradient construction depending on robust
                % optimization
                if strcmp(cst{i,6}(j).robustness,'none')
                
                    d_i = d{1}(cst{i,4}{1});

                    delta{1}(cst{i,4}{1}) = delta{1}(cst{i,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                end
                
            end
       
        end
            
    end
    
end
  
% Calculate gradient
g = zeros(dij.totalNumOfBixels,1);

for i = 1:options.numOfScenarios
    
    if any(delta{i} ~= 0) % exercise only if contributions from scenario i
        
        if isequal(options.bioOpt,'none')
            
            if isfield(dij,'optBixel')
                g(dij.optBixel) = g(dij.optBixel) + (delta{i}' * dij.physicalDose{i}(:,dij.optBixel))';
                
                if dij.memorySaver
                    depthOffset = uint32(0);
                    tailOffset = uint32(0);
                    
                    for j = 1:dij.totalNumOfRays
                        if ~dij.optBixel(j)
                            continue
                        end
                        depthInd = depthOffset+(1:uint32(dij.nDepth(j)));
                        depthOffset = depthOffset+uint32(dij.nDepth(j));
                        
                        for k = depthInd
                            tailInd = tailOffset+(1:uint32(dij.nTailPerDepth(k)));
                            tailOffset = tailOffset+uint32(dij.nTailPerDepth(k));
                            
                            voxInd = dij.ixTail(tailInd);
                            
                            g(j) = g(j)+sum(delta{i}(voxInd)).*dij.bixelDoseTail(k);
                        end
                    end
                end
                
            else
                g = g + (delta{i}' * dij.physicalDose{i})';
                
                if dij.memorySaver
                    depthOffset = uint32(0);
                    tailOffset = uint32(0);
                    
                    for j = 1:dij.totalNumOfRays
                        depthInd = depthOffset+(1:uint32(dij.nDepth(j)));
                        depthOffset = depthOffset+uint32(dij.nDepth(j));
                        
                        for k = depthInd
                            tailInd = tailOffset+(1:uint32(dij.nTailPerDepth(k)));
                            tailOffset = tailOffset+uint32(dij.nTailPerDepth(k));
                            
                            voxInd = dij.ixTail(tailInd);
                            
                            g(j) = g(j)+sum(delta{i}(voxInd)).*dij.bixelDoseTail(k);
                        end
                    end
                end
            end
            
            g = g.*dij.scaleFactor;

        elseif isequal(options.ID,'protons_const_RBExD')
            
            g            = g + (delta{i}' * dij.physicalDose{i} * dij.RBE)';
            
        elseif isequal(options.bioOpt,'LEMIV_effect')

            vBias        = (delta{i}' * dij.mAlphaDose{i})';
            quadTerm     = dij.mSqrtBetaDose{i} * w;
            mPsi         = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{i})';
            g            =  g + vBias + mPsi ; 

        elseif isequal(options.bioOpt,'LEMIV_RBExD')

            deltaTmp              = zeros(dij.numOfVoxels,1);
            scaledEffect          = d{i} + dij.gamma;
            deltaTmp(dij.ixDose)  = delta{i}(dij.ixDose)./(2*dij.bx(dij.ixDose).*scaledEffect(dij.ixDose));
            vBias                 = (deltaTmp' * dij.mAlphaDose{i})';
            quadTerm              = dij.mSqrtBetaDose{i} * w;
            mPsi                  = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{i})';
            g                     = g + vBias + mPsi ;

        end

    end
end
