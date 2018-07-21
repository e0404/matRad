function PSIrange =  matRad_calcSecRangeMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m,w_j,w_m,mW_CovBio) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the second raw moment. Please note that i
% and l depict voxel indices and j and m pencil beam indices. This function
% allows to calculate the the correlation between on spot j and multiple
% spots m simultaniously. 
% 
% call
%   PSI_ijlm =  matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m)
%
% input
%   vLaSi11:        Lambda + Sigma of spot combination j,j 
%   vLaSi22:        Lambda + Sigma of spot combination m,m  
%   vLaSi12:        Lambda + Sigma of spot combination j,m 
%   vLaSi21:        Lambda + Sigma of spot combination m,j 
%   Dev_j:          distance between radiological depth and the Gaussian means of component j
%   Dev_m:          distance between radiological depth and the Gaussian means of component m
%   mW_CovBio:      covariance of the gaussian weights
%
% output
%   PSI_ijlm:         matRad critical structure struct cst
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(isa(vLaSi11,'double'));
assert(isa(vLaSi22,'double'));  
assert(isa(vLaSi12,'double'));  
assert(isa(vLaSi21,'double'));  
assert(isa(Dev_j,'double'));  
assert(isa(Dev_m,'double'));  
assert(isa(w_j,'double'));  
assert(isa(w_m,'double'));  
assert(isa(mW_CovBio,'double'));  

upperBound = 20;

coder.varsize('vLaSi11',[upperBound 1],[1 0]);
coder.varsize('vLaSi22',[upperBound 1],[1 0]);
coder.varsize('vLaSi12',[1 1]);
coder.varsize('vLaSi21',[1 1]);
coder.varsize('Dev_j',[upperBound 1],[1 0]);
coder.varsize('Dev_m',[upperBound 1],[1 0]);
coder.varsize('w_j',[upperBound 1],[1 0]);
coder.varsize('w_m',[upperBound 1],[1 0]);
coder.varsize('mW_CovBio',[upperBound upperBound],[1 1]);


bioOffset = mW_CovBio + (w_j * w_m');
Det       = vLaSi11*vLaSi22' - (vLaSi12*vLaSi21');
FracDet   = 1./(2*pi*sqrt(Det));
ExpTerm   = FracDet .* exp(-.5./Det.*((Dev_j.^2)*vLaSi22' - bsxfun(@times,(Dev_j*Dev_m'),(vLaSi12+vLaSi21)) + vLaSi11*(Dev_m.^2)'));
PSIrange  = ones(numel(w_j),1)' * (ExpTerm .* bioOffset) * ones(numel(w_m),1);
 
end





