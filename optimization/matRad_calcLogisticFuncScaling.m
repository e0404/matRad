function k = matRad_calcLogisticFuncScaling(x,xRef,xRatio,refVal,kMin,kMax)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad calculation of parameter k in logistic function f(x) = 1/(1+exp(-2k[x-xRef]))
% 
% call
%   k = matRad_calcLogisticFuncScaling(x,xRef,xRatio,refVal,kMin,kMax)
%
% input
%   x      : x values
%   xRef   : reference x value
%   xRatio : ratio of x values that should be included in [xRef - xDelta,xRef + xDelta]
%   refVal : reference value of logistic function such that min(f(x)) = refVal
%   kMin   : minimum value of scaling parameter
%   kMax   : maximum value of scaling parameter
%
% output
%   k:  scaling parameter k in logistic function
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

% calculate absolute differences from xRef, sort absolute differences
absDiffSorted = sort(abs(xRef-x));

% determine x interval [xRef - xDelta,xRef + xDelta] such that xRatio*numel(x) lies in this interval
xDelta = absDiffSorted(round(xRatio*numel(x)));

% calculate scaling parameter k such that min(f(x)) = refVal or max(f(x)) = 1 - refVal
k = log(1/refVal-1)/(2*xDelta);

% apply minimum and maximum limits for k
k = max(k,kMin);
k = min(k,kMax);


end