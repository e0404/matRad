function [f, g, d] = matRad_IMRTBioObjFunc(w,dij,cst)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call [f, g, d] = matRad_IMRTBioObjFunc(w,dij,cst)
% to calculate the biologic objective function value f, the gradient g, and the dose
% distribution d
% f: objective function value
% g: gradient vector
% bd: biological dose vector
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) by Mark Bangert 2014
% m.bangert@dkzf.de

% Calculate biological effect
d = dij.dose*w;
a = ((dij.mAlpha.*dij.dose)*w);
a(isnan(a))=0;
b = (dij.mBeta).* (dij.dose*w.^2);
b(isnan(b))=0;


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
        
        % Minimun penalty
        rho_min = cst{i,7};
        
        % Maximum penalty
        rho_max = cst{i,6};
        
        % get dose, alpha and beta vector in current VOI
        a_i = a(cst{i,8});
        b_i = b(cst{i,8});
        
        % Maximun deviation: biologic effect minus maximun prescribed biological effect.
        deviation_max = a_i + b_i - cst{i,4};
        
        % Minimun deviation: Dose minus minimun dose.
        deviation_min = a_i + b_i - cst{i,5};
        
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

% Calculate gradient.
%g = 2 * (delta' * dij.dose)';
a_w = dij.mAlpha*w;
b_w = 2.*(dij.mBeta.* dij.dose*w);
b_w(isnan(b_w))=0;

vTmp = a_w+b_w;
g = 2 * ((delta.*vTmp)'*dij.dose)';

%g = 2 * ((delta.*(a + 2*b.*d))'*dij.dose)';

%% frist gradient from weight w1  
vFac=(dij.mAlpha(:,1)+0.1.*dij.dose*w);
g1 = (delta.*vFac)' *dij.dose(:,1);

%% calculate all gradients
vec = 0.1.*dij.dose*w;
n = length(0.1.*dij.dose*w);
mDia= spdiags(vec(:),0,n,n);
mColumnMulti = mDia*dij.mAlpha;
mInnerDeviation = mColumnMulti.*dij.dose;

g2=  2*(delta'*mInnerDeviation)';


M = sparse(vec * (ones(size(dij.mAlpha,2),1))');





