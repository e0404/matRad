function vY =  matRad_sumGauss(vX,vMu,vSqSigma,vW) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the sum of Gaussians
% 
% call
%   vY =  matRad_sumGauss(vX,vMu,vSqSigma)
%
% input
%   vX:         query points at which the sum of gaussians should be
%               evaluated
%   vMu:        mean vector
%   vSqSigma:   squared sigma vector
%   vW          weights
%
% output
%   vY:         sum of gaussians at vX
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

assert(isa(vX,'double'));
assert(isa(vMu,'double'));  
assert(isa(vSqSigma,'double'));  
assert(isa(vW,'double')); 

NumComp = 50;
coder.varsize('vX',[Inf 1],[1 0]);
coder.varsize('vMu',[NumComp 1],[1 0]);
coder.varsize('vSqSigma',[NumComp 1],[1 0]);
coder.varsize('vW',[NumComp 1],[1 0]);

vY =  ( exp(-bsxfun(@minus,vX,vMu').^2 ./ (2* ones(numel(vX),1) * vSqSigma' ))./(sqrt(2*pi*ones(numel(vX),1) * vSqSigma'))) * vW ;
                                
end

