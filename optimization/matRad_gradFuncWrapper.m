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

% get current dose / effect / RBExDose vector
[d,d_exp,Omega] = matRad_backProjection(w,dij,cst,options);

% Initializes delta
delta         = cell(options.numOfScen,1);
[delta{:}]    = deal(zeros(dij.numOfVoxels,1));
delta_exp{1}  = zeros(dij.numOfVoxels,1);
vOmega        = 0;

% if composite worst case optimization is used then create a cell array for book keeping
for i = 1:size(cst,1)
  for j = 1:numel(cst{i,6})
      if strcmp(cst{i,6}(j).robustness,'COWC')
         f_COWC          = zeros(options.numOfScen,1);
         delta_COWC      = cell(options.numOfScen,1);
        [delta_COWC{:}]  = deal(zeros(dij.numOfVoxels,1)); break;
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
                if (~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')) &&...
                    isequal(options.quantityOpt,'effect') 

                    d_ref = cst{i,5}.alphaX*cst{i,6}(j).dose + cst{i,5}.betaX*cst{i,6}(j).dose^2;
                else
                    d_ref = cst{i,6}(j).dose;
                end
                
                % different gradient construction depending on robust
                % optimization
                if strcmp(cst{i,6}(j).robustness,'none')
                
                    d_i = d{1}(cst{i,4}{1});

                    delta{1}(cst{i,4}{1}) = delta{1}(cst{i,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);
                    
                elseif strcmp(cst{i,6}(j).robustness,'PROB')

                        d_i = d_exp{1}(cst{i,4}{1});

                        delta_exp{1}(cst{i,4}{1}) = delta_exp{1}(cst{i,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);

                        vOmega = vOmega + Omega{i};
                        
                 % VWWC = voxel-wise worst case; VWWC_CONF = conformity voxel-wise worst case    
                 elseif strcmp(cst{i,6}(j).robustness,'VWWC') || strcmp(cst{i,6}(j).robustness,'VWWC_CONF')

                    % prepare min/max dose vector for voxel-wise worst case  
                    if ~exist('d_tmp','var')
                         d_tmp = [d{:}];
                    end   
                    
                    d_Scen = d_tmp(cst{i,4}{1},:);
                    [d_max,max_ix] = max(d_Scen,[],2);
                    [d_min,min_ix] = min(d_Scen,[],2);
                    
                    if isequal(cst{i,3},'OAR')
                        d_i = d_max;
                    elseif isequal(cst{i,3},'TARGET')          
                        d_i = d_min;
                    end

                    if sum(isnan(d_min)) > 0
                        warning('nan values in gradFuncWrapper');
                    end
                    
                    if strcmp(cst{i,6}(j).robustness,'VWWC')
                        deltaTmp = matRad_gradFunc(d_i,cst{i,6}(j),d_ref);
                    elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square overdosing')
                        deltaTmp = matRad_gradFunc(d_max,cst{i,6}(j),d_ref);
                    elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square underdosing')
                        deltaTmp = matRad_gradFunc(d_min,cst{i,6}(j),d_ref);
                       
                    end
                    
                    for ixScen = 1:options.numOfScen

                       if strcmp(cst{i,6}(j).robustness,'VWWC')
                           
                           if isequal(cst{i,3},'OAR')
                               currWcIx = max_ix == ixScen;   
                           elseif isequal(cst{i,3},'TARGET')
                               currWcIx = min_ix == ixScen;
                           end
                        
                           delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + deltaTmp.*currWcIx;
                           
                        elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square overdosing')
                          
                               currWcIx = max_ix == ixScen;
                               delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + deltaTmp.*currWcIx;
                           
                        elseif strcmp(cst{i,6}(j).robustness,'VWWC_CONF') && isequal(cst{i,6}(j).type, 'square underdosing')
                          
                               currWcIx = min_ix == ixScen;
                               delta{ixScen}(cst{i,4}{1}) = delta{ixScen}(cst{i,4}{1}) + deltaTmp.*currWcIx;
                           
                        end

                    end 
                    
                 % composite worst case consideres ovarall the worst objective function value   
                 elseif strcmp(cst{i,6}(j).robustness,'COWC')
                   
                       for ixScen = 1:options.numOfScen

                           d_i = d{ixScen}(cst{i,4}{1});

                           f_COWC(ixScen)                     = f_COWC(ixScen) + matRad_objFunc(d_i,cst{i,6}(j),d_ref);
                           delta_COWC{ixScen}(cst{i,4}{1})    = delta_COWC{ixScen}(cst{i,4}{1}) + matRad_gradFunc(d_i,cst{i,6}(j),d_ref);
                       end   
                       
                 % objective-wise worst case consideres the worst individual objective function value        
                 elseif strcmp(cst{i,6}(j).robustness,'OWC')
                     
                     matRad_dispToConsole(['not yet implemented \n'],param,'error');
                     
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

for i = 1:options.numOfScen
    if any(delta{i} ~= 0) || any(delta_exp{1} ~= 0) % exercise only if contributions from scenario i

        if isequal(options.quantityOpt,'physicalDose')

            g            = g + (delta{i}' * dij.physicalDose{options.ixForOpt(i)})';
            if i == 1
                g        = g + ((delta_exp{i}' * dij.physicalDoseExp{1})' + (2 * vOmega)); 
            end

        elseif isequal(options.model,'constRBE')
            
            g            = g + (delta{i}' * dij.physicalDose{options.ixForOpt(i)} * dij.RBE)';
            if i == 1
                g        = g + (delta_exp{i}' * dij.physicalDoseExp{options.ixForOpt(i)} * dij.RBE)' + (2 * vOmega);
            end
            
        elseif isequal(options.quantityOpt,'effect')

            vBias        = (delta{i}' * dij.mAlphaDose{options.ixForOpt(i)})';
            quadTerm     = dij.mSqrtBetaDose{options.ixForOpt(i)} * w;
            mPsi         = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{options.ixForOpt(i)})';
            g            =  g + vBias + mPsi ; 
            
            if i == 1
                vBias        = (delta_exp{i}' * dij.mAlphaDoseExp{1})';
                quadTerm     = dij.mSqrtBetaDoseExp{1} * w;
                mPsi         = (2*(delta_exp{i}.*quadTerm)'*dij.mSqrtBetaDoseExp{1})';
                g            =  g + vBias + mPsi ; 
            end

        elseif isequal(options.quantityOpt,'RBExD')

            deltaTmp              = zeros(dij.numOfVoxels,1);
            scaledEffect          = d{i} + dij.gamma;
            deltaTmp(dij.ixDose)  = delta{i}(dij.ixDose)./(2*dij.betaX(dij.ixDose).*scaledEffect(dij.ixDose));
            vBias                 = (deltaTmp' * dij.mAlphaDose{options.ixForOpt(i)})';
            quadTerm              = dij.mSqrtBetaDose{options.ixForOpt(i)} * w;
            mPsi                  = (2*(delta{i}.*quadTerm)'*dij.mSqrtBetaDose{options.ixForOpt(i)})';
            g                     = g + vBias + mPsi ;
            
            if i == 1
                deltaTmp              = zeros(dij.numOfVoxels,1);
                scaledEffect          = d_exp{1} + dij.gamma;
                deltaTmp(dij.ixDose)  = delta_exp{i}(dij.ixDose)./(2*dij.betaX(dij.ixDose).*scaledEffect(dij.ixDose));
                vBias                 = (deltaTmp' * dij.mAlphaDoseExp{1})';
                quadTerm              = dij.mSqrtBetaDoseExp{1} * w;
                mPsi                  = (2*(delta_exp{i}.*quadTerm)'*dij.mSqrtBetaDoseExp{1})';
                g                     = g + vBias + mPsi ;
            end
        end

    end
end
