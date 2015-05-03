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
delta = zeros(numVoxels,1);

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
                delta(cst{i,4}) = delta(cst{i,4}) + (rho/size(cst{i,4},1))*underdose;
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                % overdose : Dose minus prefered dose
                overdose = e_i - e_ref;
                
                % apply positive operator
                overdose(overdose<0) = 0;
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*(overdose'*overdose);
                
                %calculate delta
                delta(cst{i,4}) = delta(cst{i,4}) + (rho/size(cst{i,4},1))*overdose;
                
            elseif isequal(cst{i,6}(j).type, 'square deviation')
                
                % deviation : Dose minus prefered dose
                deviation = e_i - e_ref;
                
                % claculate objective function
                f = f + (rho/size(cst{i,4},1))*(deviation'*deviation);
                
                % calculate delta
                delta(cst{i,4}) = delta(cst{i,4}) + (rho/size(cst{i,4},1))*deviation;
                
            else
                
                error('undefined objective in cst struct');
                
            end

        end
        
    end
end

% gradient calculation
if nargout > 1
    vBias= (delta' * dij.mAlphaDose)';
    mPsi = ((delta.*quadTerm)'*dij.physicalDose)';
    g = 2*(vBias+mPsi);    
    
    % first gradient 
    vTemp = dij.mAlphaDose(:,1)+(dij.mSqrtBetaDose*w);
    g1 = (delta.*vTemp)'*dij.physicalDose(:,1);
     fprintf(['first gradient ' num2str(g1) ]);
end
