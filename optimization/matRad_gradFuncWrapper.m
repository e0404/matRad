function g = matRad_gradFuncWrapper(w,dij,cst,type)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: gradient function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage,
% squared deviation and DVH objectives
% 
% call
%   g = matRad_gradFuncWrapper(w,dij,cst,type)
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

global fScaling

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,type);

% Initializes delta
delta      = cell(dij.numOfScenarios,1);
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
                    isequal(type,'effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % different gradient construction depending on robust
                % optimization
                if strcmp(cst{i,6}(j).robustness,'none')
                
                    d_i = d{1}(cst{i,4}{1});

                    delta{1}(cst{i,4}{1}) = delta{1}(cst{i,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                elseif strcmp(cst{i,6}(j).robustness,'probabilistic')

                    for k = 1:dij.numOfScenarios

                        d_i = d{k}(cst{i,4}{1});

                        delta{k}(cst{i,4}{1}) = delta{k}(cst{i,4}{1}) + dij.probOfScenarios(k)*matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                    end

                elseif strcmp(cst{i,6}(j).robustness,'voxel-wise worst case')

                    % prepare min/max dose vector we have chosen voxel-wise worst case
                    if ~exist('d_max','var')
                        [d_max,max_ix] = max([d{:}],[],2);
                        [d_min,min_ix] = min([d{:}],[],2);
                    end

                    if isequal(cst{i,3},'OAR')
                        d_i = d_max(cst{i,4}{1});
                    elseif isequal(cst{i,3},'TARGET')
                        d_i = d_min(cst{i,4}{1});
                    end

                    deltaTmp = matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                    for k = 1:dij.numOfScenarios

                        if isequal(cst{i,3},'OAR')
                            currWcIx = max_ix(cst{i,4}{1}) == k;

                        elseif isequal(cst{i,3},'TARGET')
                            currWcIx = min_ix(cst{i,4}{1}) == k;
                        end

                        delta{k}(cst{i,4}{1}) = delta{k}(cst{i,4}{1}) + deltaTmp.*currWcIx;

                    end
                    
                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    % calc invers DCH
                    Q_ref  = cst{i,6}(j).coverage/100;
                    V_ref  = cst{i,6}(j).volume/100;
                    d_ref2 = matRad_calcInversDCH(V_ref,Q_ref,d,dij,cst(i,:)); 
                    
                    if dij.numOfScenarios > 1
                        
                        for k = 1:dij.numOfScenarios
                            
                            % get VOI dose in current scenario
                            d_i = d{k}(cst{i,4}{1});

                            % get voxel dependent weigthing
                            voxelWeighting = 1; 
                            
                            % calculate dose deviations from d_ref
                            delta{k}(cst{i,4}{1}) = delta{k}(cst{i,4}{1}) + dij.ScenProb(k)*matRad_gradFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);    
                                                                 
                        end
                        
                    else
                        
                        % get VOI ScenUnion dose of nominal scneario
                        cstLogical = strcmp(cst(:,2),[cst{i,2},' ScenUnion']);
                        d_i        = d{1}(cst{cstLogical,5}.voxelID);
                        
                        % get voxel dependent weigthing
                        voxelWeighting = 5*cst{cstLogical,5}.voxelProb;  
                        
                        % calculate delta
                        delta{1}(cst{cstLogical,4}{1}) = delta{1}(cst{cstLogical,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);
                        
                    end
                    
                end
                
            end
       
        end
            
    end
    
end

  
% Calculate gradient
g = zeros(dij.totalNumOfBixels,1);

for i = 1:dij.numOfScenarios
    if any(delta{i} > 0) % exercise only if contributions from scenario i

        if isequal(type,'none')

            g = g + (delta{i}' * dij.physicalDose{i})';

        elseif isequal(type,'effect')

            vBias    = (delta{i}' * dij.mAlphaDose{i})';
            quadTerm = dij.mSqrtBetaDose{i} * w;
            mPsi     = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{i})';
            g        =  g + vBias + mPsi ; 

        elseif isequal(type,'RBExD')

            scaledEffect = d{i} + dij.gamma;
            deltaTmp     = delta{i}./(2*dij.bx.*scaledEffect);
            vBias        = (deltaTmp' * dij.mAlphaDose{i})';
            quadTerm     = dij.mSqrtBetaDose{i} * w;
            mPsi         = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{i})';
            g            = g + vBias + mPsi ;

        end

    end
end

% apply objective scaling
g = fScaling.*g;

% save min/max gradient
global matRad_iteration
global GRADIENT
GRADIENT(1,1,matRad_iteration+1)= max(abs(g));
GRADIENT(1,2,matRad_iteration+1)= min(abs(g));