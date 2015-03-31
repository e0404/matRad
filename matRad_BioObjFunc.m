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

% default alpha photon and beta photon parameters to calculate prescribed effect
a_x = 0.1; 
b_x = 0.05;

% Numbers of voxels
numVoxels = size(dij.dose,1);

% Initializes f
f = 0;

% Initializes delta
delta = zeros(numVoxels,1);

% Compute optimization function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
        
        % get tissue specific alpha photon and beta photon to calculate
        % prescriped effect
        if ~isempty(cst{i,9})
            a_x = cst{i,9}.alphaX;
            b_x = cst{i,9}.betaX;
        end
        
        % prescribed effect
        Emax = a_x*cst{i,4}+b_x*cst{i,4}^2;
        Emin = a_x*cst{i,5}+b_x*cst{i,5}^2;

        % Minimun penalty
        rho_min = cst{i,7};
        
        % Maximum penalty
        rho_max = cst{i,6};
        
        % get effect vector in current VOI
        e_i = e(cst{i,8});
        
        % Maximun deviation: biologic effect minus maximun prescribed biological effect.
        deviation_max = e_i - Emax;
        
        % Minimun deviation: Dose minus minimun dose.
        deviation_min = e_i - Emin;
        
        % Apply positive operator H.
        deviation_max(deviation_max<0) = 0;
        deviation_min(deviation_min>0) = 0;
                
        % Calculate the objective function
        f = f + (rho_max/size(cst{i,8},1))*(deviation_max'*deviation_max) + ...
            (rho_min/size(cst{i,8},1))*(deviation_min'*deviation_min);
        
        % Calculate delta
        delta(cst{i,8}) = delta(cst{i,8}) + (rho_max/size(cst{i,8},1))*deviation_max +...
            (rho_min/size(cst{i,8},1))*deviation_min;

    end
end

% gradient calculation - compromises one expensive calculation
if nargout > 1
    profile on
   
    vBias= (delta' * dij.mAlphaDose)';
    %mPsi= (delta'*((spdiags(2*quadTerm(:),0,n,n)*dij.doseSkeleton).*dij.dose))';
    mPsi = ((2*delta.*quadTerm)'*dij.dose)';
    g = 2*(vBias+mPsi);    
end
