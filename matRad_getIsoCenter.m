function isoCenter = matRad_getIsoCenter(cst,ct,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the isocenter [mm] as the joint center of gravity of all volumes
% of interest that are labeled as target within the cst struct
% 
% call
%   isoCenter = matRad_getIsoCenter(cst,ct,pln,visBool)
%
% input
%   cst:        matRad cst struct
%   ct:         ct cube
%   visBool:    toggle on/off visualization (optional)
%
% output
%   isoCenter:  isocenter in [mm]   
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

% if visBool not set toogle off visualization
if nargin < 3
    visBool = 0;
end

% Initializes V variable.
V = [];

% Save target indices in V variable.
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;cst{i,4}];
    end
end

% Delete repeated indices, one voxel can belong to two VOIs, because
% VOIs can be overlapping.
V = unique(V);

% throw error message if no target is found
if isempty(V)
    error('Could not find target');
end

% Transform subcripts from linear indices 
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct.cube),V);

% Transform to [mm]
xCoordsV = xCoordsV * ct.resolution(1);
yCoordsV = yCoordsV * ct.resolution(2);
zCoordsV = zCoordsV * ct.resolution(3);

% Calculated isocenter.
isoCenter = mean([xCoordsV yCoordsV zCoordsV]);

% Visualization
if visBool

    clf
    hold on
    
    % Plot target
    plot3(yCoordsV,xCoordsV,zCoordsV,'kx')
    
    % Show isocenter: red point
    plot3(isoCenter(2),isoCenter(1),isoCenter(3),'r.','MarkerSize',30)
    
    xlabel('y [mm]')
    ylabel('x [mm]')
    zlabel('z [mm]')
    
end
