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

global fScaling

% get current dose / effect / RBExDose vector
d = matRad_backProjection(w,dij,options);

% Initializes delta
delta      = cell(options.numOfScenarios,1);
[delta{:}] = deal(zeros(dij.numOfVoxels,1));

for i = 1:size(cst,1)
  for j = 1:numel(cst{i,6})
      if strcmp(cst{i,6}(j).robustness,'COWC')
         f_COWC = zeros(options.numOfScenarios,1);
         delta_COWC      = cell(options.numOfScenarios,1);
        [delta_COWC{:}]  = deal(zeros(dij.numOfVoxels,1));
      end
  end
end

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints and objectives for the current VOI
        for j = 1:numel(cst{i,6})
        
            % only perform gradient computations for objectives
            if isempty(strfind(cst{i,6}(j).type,'constraint'))

                % compute reference
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) && isequal(options.quantity,'effect')

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

                    for ixScen = 1:dij.numOfScenarios

                        d_i = d{ixScen}(cst{i,4}{1});

                        delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + dij.probOfScenarios(ixScen) * matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                    end

                elseif strcmp(cst{i,6}(j).robustness,'VWWC')

                    % prepare min/max dose vector for voxel-wise worst case
                    if ~exist('d_max','var')
                        [d_max,max_ix] = max([d{:}],[],2);   
                        [d_min,min_ix] = min([d{:}],[],2);
                    end

                    if isequal(cst{i,3},'OAR')
                        d_i = d_max(cst{i,4}{1});
                    elseif isequal(cst{i,3},'TARGET')
                        d_i = d_min(cst{i,4}{1});
                    end

                    if sum(isnan(d_min)) > 0
                        warning('nan values in gradFuncWrapper');
                    end
                    
                    deltaTmp = matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                    for ixScen = 1:dij.numOfScenarios

                        if isequal(cst{i,3},'OAR')
                            currWcIx = max_ix(cst{i,4}{1}) == ixScen;   %

                        elseif isequal(cst{i,3},'TARGET')
                            currWcIx = min_ix(cst{i,4}{1}) == ixScen;
                        end

                        delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + deltaTmp.*currWcIx;

                    end
                 
                elseif strcmp(cst{i,6}(j).robustness,'COWC')
                   
                    for ixScen = 1:dij.numOfScenarios

                        d_i = d{ixScen}(cst{i,4}{1});

                        f_COWC(ixScen) = f_COWC(ixScen) + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                        delta_COWC{ixScen}(cst{i,4}{1}) = delta_COWC{ixScen}(cst{i,4}{1}) +  matRad_gradFunc(d_i,cst{i,6}(j),d_ref);
                    end
                     

                elseif strcmp(cst{i,6}(j).robustness,'coverage')
                    
                    % calc invers DCH
                    Q_ref  = cst{i,6}(j).coverage/100;
                    V_ref  = cst{i,6}(j).volume/100;
                    d_ref2 = matRad_calcInversDCH(V_ref,Q_ref,d,dij,cst(i,:)); 
                    
                    if dij.numOfScenarios > 1
                        
                        for ixScen = 1:dij.numOfScenarios
                            
                            % get VOI dose in current scenario
                            d_i = d{ixScen}(cst{i,4}{1});

                            % get voxel dependent weigthing
                            voxelWeighting = 1; 
                            
                            % calculate dose deviations from d_ref
                            delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + dij.ScenProb(ixScen)*matRad_gradFunc(d_i,cst{i,6}(j),d_ref,d_ref2,voxelWeighting);    
                                                                 
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


% extract current worst case scenario
if exist('f_COWC','var')
   [~,ixCurrWC]    = max(f_COWC(:));
   delta{ixCurrWC} = delta_COWC{ixCurrWC};
end

% Calculate gradient
g = zeros(dij.totalNumOfBixels,1);

for i = 1:options.numOfScenarios
    if any(delta{i}) % exercise only if contributions from scenario i

        if isequal(options.quantity,'physicalDose')

            g            = g + (delta{i}' * dij.physicalDose{dij.indexforOpt(i)})';

        elseif isequal(options.type,'const_RBExD')
            
            g            = g + (delta{i}' * dij.physicalDose{dij.indexforOpt(i)} * dij.RBE)';
            
        elseif isequal(options.quantity,'effect') 

            vBias        = (delta{i}' * dij.mAlphaDose{dij.indexforOpt(i)})';
            quadTerm     = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
            mPsi         = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{dij.indexforOpt(i)})';
            g            =  g + vBias + mPsi ; 

        elseif isequal(options.quantity,'RBExD') 

            scaledEffect = d{i} + dij.gamma;
            deltaTmp     = delta{i}./(2*dij.bx.*scaledEffect);
            vBias        = (deltaTmp' * dij.mAlphaDose{dij.indexforOpt(i)})';
            quadTerm     = dij.mSqrtBetaDose{dij.indexforOpt(i)} * w;
            mPsi         = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{dij.indexforOpt(i)})';
            g            = g + vBias + mPsi ;
            
            % utils
            ab = dij.ax./dij.bx;
            sqab = real(sqrt(ab));
            dp   = dij.physicalDose{dij.indexforOpt(i)} * w;
            LETd = (dij.mLETDose{dij.indexforOpt(i)} * w)./dp;
            
            % deriviative effect
%             NumGrad = 1;
%             A =  options.p0 .* dij.ax .* dij.physicalDose{1}(:,NumGrad) + options.p1 * dij.bx .* dij.mLETDose{1}(:,NumGrad);
%             B = 2 * dij.bx .* (options.p2 * dp - options.p3 * sqab .* LETd .* dp) .* (options.p2 * dij.physicalDose{1}(:,NumGrad) - options.p3 * sqab .* dij.mLETDose{1}(:,NumGrad));
%             
%             gnew = delta{i}'*(A + B);
            
%             
%             % 
%             RBEmax = options.p0 + ((options.p1 ./ab).*LETd);
%             RBEmin = options.p2 - (options.p3 * sqrt(ab).*LETd);
%              
%             Fac   = 1 ./(sqrt(ab.^2 + (4*dp.*ab.*RBEmax) + (4*dp.^2 .* RBEmin.^2))) ;
%             
%             for NumGrad = 1:7
% 
%                part1 = (4.*ab.* options.p0.*dij.physicalDose{1}(:,NumGrad)) + (4 * options.p1 * dij.mLETDose{1}(:,NumGrad));
% 
%                part2 = 8*(options.p2*dp + options.p3 * sqab.* (dij.mLETDose{dij.indexforOpt(i)} * w)).* (options.p2.*dij.physicalDose{1}(:,NumGrad) + options.p3 .* sqab .* dij.mLETDose{1}(:,NumGrad));
% 
%                A = Fac .* (part1 + part2);
%                A(isnan(A)) = 0;
% 
%                g_single = 0.5*(A'*delta{i});
% 
%                g(NumGrad) = g_single;
%             
%             end

            
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