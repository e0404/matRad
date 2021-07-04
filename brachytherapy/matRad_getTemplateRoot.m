function templateRoot = matRad_getTemplateRoot(ct,cst)
%MATRADGETTEMPLATEROOT calculates origin position for template
% INPUT
%   ct:         ct cube
%   cst:        matRad cst struct
% OUTPUT
%   templateRoot: 
%               1x3 column vector with root position
%               x,y : center \\ z : bottom of target VOI

% if visBool not set toogle off visualization
if nargin < 3
    visBool = 0;
end
%% Initializes V variable.
V = [];

%Check if any constraints/Objectives have been defined yet
noObjOrConst = all(cellfun(@isempty,cst(:,6)));

% Save target indices in V variable.
for i = 1:size(cst,1)
    % We only let a target contribute if it has an objective/constraint or
    % if we do not have specified objectives/constraints at all so far
    if isequal(cst{i,3},'TARGET') && (~isempty(cst{i,6}) || noObjOrConst)
        V = [V; cst{i,4}{1}];
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
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(ct.cubeDim,V);

%% Transform to [mm] in ct grid
xCoordsV = ct.x(xCoordsV);
yCoordsV = ct.y(yCoordsV);
zCoordsV = ct.z(zCoordsV);

%% define root position
templateRoot = [mean(xCoordsV),mean(yCoordsV),min(zCoordsV)];

end

