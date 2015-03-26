function [f, g, d, bd] = matRad_IMRTObjFunc(w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad objective function for inverse planning
% 
% call
%   [f, g, d] = matRad_IMRTObjFunc(w,dij,cst)
%
% input
%   w:   bixel weight vector
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   f: objective function value
%   g: gradient
%   d: dose distribution
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

bd=NaN;

% Calculate dose
d = dij*w;

% Numbers of voxels
numVoxels = size(dij,1);

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
        
        % get dose vector in current VOI
        d_i = d(cst{i,8});
        
        % Maximun deviation: Dose minus maximun dose.
        deviation_max = d_i - cst{i,4};
        
        % Minimun deviation: Dose minus minimun dose.
        deviation_min = d_i - cst{i,5};
        
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
if nargout > 1
    g = 2 * (delta' * dij)';
end
