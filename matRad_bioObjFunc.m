function [f, g] = matRad_bioObjFunc(w,dij,cst)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call [f, g, d] = matRad_IMRTBioObjFunc(w,dij,cst)
% to calculate the biologic objective function value f, the gradient g, and the dose
% distribution d
% f: objective function value
% g: gradient vector
% bd: biological effect vector
% d: physical dose vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) by Mark Bangert 2014
% m.bangert@dkzf.de



% calculate biological effect
linTerm  = dij.mAlphaDose*w;
quadTerm = dij.mSqrtBetaDose*w;
e = linTerm + quadTerm.^2;

% Numbers of voxels
numVoxels = size(dij.physicalDose,1);

% Initializes f
f = 0;

% Initializes delta
delta_underdose = zeros(numVoxels,1);
delta_overdose  = zeros(numVoxels,1);
delta_deviation = zeros(numVoxels,1);
delta_mean      = zeros(numVoxels,1);
delta_EUD       = zeros(numVoxels,1);



% Compute optimization function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
        
        % get effect vector in current VOI
        e_i = e(cst{i,4});
        
        % get tissue specific alpha photon and beta photon to calculate
        % prescriped effect
        a_x = cst{i,5}.alphaX;
        b_x = cst{i,5}.betaX;
        
        % loop over the number of constraints for the current VOI
        for j = 1:size(cst{i,6},2)
            
            % get Penalty
            rho = cst{i,6}(j).parameter(1);
            
            % refernce effect
            e_ref = a_x*cst{i,6}(j).parameter(2)+b_x*cst{i,6}(j).parameter(2)^2;
            
            if isequal(cst{i,6}(j).type, 'square underdosing')
  
                % underdose : effect minus reference effect
                underdose = e_i - e_ref;

                % apply positive operator
                underdose(underdose>0) = 0;
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*(underdose'*underdose);
                
                % calculate delta
                delta_underdose(cst{i,4}) = delta_underdose(cst{i,4}) + (rho/size(cst{i,4},1))*underdose;
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                % overdose : Dose minus prefered dose
                overdose = e_i - e_ref;
                
                % apply positive operator
                overdose(overdose<0) = 0;
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*(overdose'*overdose);
                
                %calculate delta
                delta_overdose(cst{i,4}) = delta_overdose(cst{i,4}) + (rho/size(cst{i,4},1))*overdose;
                
            elseif isequal(cst{i,6}(j).type, 'square deviation')
                
                % deviation : Dose minus prefered dose
                deviation = e_i - e_ref;
                
                % claculate objective function
                f = f + (rho/size(cst{i,4},1))*(deviation'*deviation);
                
                % calculate delta
                delta_deviation(cst{i,4}) = delta_deviation(cst{i,4}) + (rho/size(cst{i,4},1))*deviation;
            
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*sum(e_i);
                
                % calculate delta
                delta_mean(cst{i,4}) = delta_mean(cst{i,4}) + ...
                    (rho/size(cst{i,4},1))*ones(size(cst{i,4},1),1);
                
             elseif isequal(cst{i,6}(j).type, 'EUD') 
                
                % get exponent for EUD
                exponent = cst{i,6}(j).exponent;
                
                % calculate objective function and delta
                if sum(e_i.^exponent)>0
                    
                    f = f + rho*nthroot((1/size(cst{i,4},1))*sum(e_i.^exponent),exponent);
                    
                    delta_EUD(cst{i,4}) = delta_EUD(cst{i,4}) + ...
                        rho*nthroot(1/size(cst{i,4},1),exponent) * sum(e_i.^exponent)^((1-exponent)/exponent) * (e_i.^(exponent-1));                    
                end    
                
                
            else
                
                error('undefined objective in cst struct');
                
            end

        end
        
    end
end

% gradient calculation
if nargout > 1
    delta = delta_underdose + delta_overdose + delta_deviation + ...
                delta_mean + delta_EUD;
    vBias= (delta' * dij.mAlphaDose)';
    mPsi = ((delta.*quadTerm)'*dij.physicalDose)';
    g = 2*(vBias+mPsi);    
    
    % first gradient - according to Paper of Jan Wilkens - Fast Multifield
%     vTemp = dij.mAlphaDose(:,1)+(dij.mSqrtBetaDose*w);
%     g1 = (2*delta.*vTemp)'*dij.physicalDose(:,1);
%     fprintf(['first gradient ' num2str(g1) '\n']);
end
