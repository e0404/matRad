function [f, g, bd] = matRad_IMRTBioObjFunc(w,dij,cst)
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

%profile on
% Calculate biological effect
d = dij.dose*w;
a = (dij.mAlphaDose*w);
%a(isnan(a))=0;
b = (dij.mBeta).* (dij.dose*w.^2);
%b(isnan(b))=0;

%biological dose
bd = a+b;

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



%% frist gradient from weight w1  
% tic
% vFac=(dij.mAlpha(:,1)+0.1.*dij.dose*w);
% g1 = 2*(delta.*vFac)' *dij.dose(:,1);
% toc

%% calculate all gradients
% vec = exp(2*dij.mBeta.*dij.dose*w);
% n = length(vec);
% %takes 0.370s
% mDia= spdiags(vec(:),0,n,n);
% % both spfun take 1.7s because two sparse matrixes are created
% mColumnExp=mDia*spfun(@exp,dij.mAlpha);
% mColumn = spfun(@log,mColumnExp);
% mInnerDev= mColumn.*dij.dose;
% g=  2*(delta'*mInnerDev)';


% if nargout > 1
%     tic
%     vVec = (2*dij.mBeta.*d);
%     n = length(vVec);
%     mA = spdiags(vVec(:),0,n,n)*dij.doseSkeleton;
%     mA = mA+dij.mAlpha;
%     mA= mA.*dij.dose;
%     g = 2*(delta'*mA)';
%     toc
% end

if nargout > 1

    lambda = (2*dij.mBeta.*d);
    n = length(lambda);
    w= (delta' * dij.mAlphaDose)';
    u= (delta'*spdiags(lambda(:),0,n,n)*dij.doseSkeleton)';  
    v= (delta'* dij.dose)';
    
    V =u.*v/sum(delta(:));
    
    g2 = 2*(V+w);

    %diff = abs(g(1)-g2(1));
end



%% repmat and bsxfun cannot be used 
%% for loop over all rows takes 0.25s per row -> *7000 ~= 30min


%% example to illustrate the idea of adding columnwise vector B
% 
% A=[4 3; 2 6];
% B = [2 3];
% 
% expA = exp(A);
% expB = exp(B);
% diaB = diag(expB);
% 
% C=diaB * expA;
% res=log(C);








