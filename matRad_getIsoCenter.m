function isoCenter = matRad_getIsoCenter(cst,ct,pln,visBool)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns the iso center of all target volumes on voi
% isoCenter: iso center [mm]

% if visBool not set toogle off visualization
if nargin < 3
    visBool = 0;
end

% Initializes V variable.
V = [];

% Save target indices in V variable.
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET')
        V = [V;cst{i,8}];
    end
end

% Delete repeated indices, one voxel can belong to two VOIs, because
% VOIs can be overlapping.
V = unique(V);

% Transform subcripts from linear indices 
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct),V);

% Transform to [mm]
xCoordsV = xCoordsV * pln.resolution(1);
yCoordsV = yCoordsV * pln.resolution(2);
zCoordsV = zCoordsV * pln.resolution(3);

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
    
    axis equal on tight
end