function  matRad_verifyGradient(func,NumBixel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to check analytically calculated gradients against finite
% differences
% 
% call
%   matRad_verifyGradient(func,NumBixel)
%
% input
%   func:     objective function handle
%   NumBixel: number of degrees of freedom
%
% output
%   The function will print the relative differences to the command line
%   prompt for 100 random components
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wInit = rand(NumBixel,1);

[f, g] = func(wInit);
epsilon = 1e-05;
NumRealizations = 100;

for i = 1:NumRealizations
    RandomComp = ceil(NumBixel*rand);
    wDelta = wInit;
    wDelta(RandomComp) = wDelta(RandomComp) + epsilon;
    [fDelta, ~] = func(wDelta);
    numGrad = (fDelta-f)/epsilon;
    diff = ((g(RandomComp)/(numGrad))-1)*100;
    fprintf(['Component # ' num2str(RandomComp) ' - percent diff in numerical and analytical gradient = '...
        num2str(diff) '\n']);
end



end

