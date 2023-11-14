function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
% matRad_DoseEngine.calcDoseInit: Interface for dose calculation
%   method for setting and preparing the inition parameters for the
%   dose calculation.
%   Should be called at the beginning of calcDose method.
%   Can be expanded or changed by overwriting this method and calling
%   the superclass method inside of it
%
% call:
%   [ct,stf,pln,dij] = matRad_DoseEngine.calcDoseInit(this,ct,stf,pln,cst)
%
% input:
%   ct:             matRad ct  struct
%   cst:            matRad cst struct
%   stf:            matRad stf struct
%
% returns:
%   dij:            matRad dij struct
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg =  MatRad_Config.instance();

this.timers.full = tic;

if this.calcDoseDirect
    msg = sprintf('Forward dose calculation using  ''%s'' Dose Engine...',this.name);
else
    msg = sprintf('Dose influence matrix calculation using  ''%s'' Dose Engine...',this.name);
end
matRad_cfg.dispInfo('%s\n',msg);

% initialize waitbar
% TODO: This should be managed from the user interface instead
if ~matRad_cfg.disableGUI
    this.hWaitbar = waitbar(0,msg);
    % prevent closure of waitbar and show busy state
    set(this.hWaitbar,'pointer','watch');
end
this.lastProgressUpdate = tic;

if numel(unique({stf(:).machine})) ~= 1 || numel(unique({stf(:).radiationMode})) ~= 1
    matRad_cfg.dispError('machine and radiation mode need to be unique within supplied stf!');
end

dij = struct();

%store CT grid
dij.ctGrid.resolution = ct.resolution;
%dij.ctGrid.resolution.x = ct.resolution.x;
%dij.ctGrid.resolution.y = ct.resolution.y;
%dij.ctGrid.resolution.z = ct.resolution.z;

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    dij.ctGrid.x = ct.resolution.x*[0:ct.cubeDim(2)-1]-ct.resolution.x/2;
    dij.ctGrid.y = ct.resolution.y*[0:ct.cubeDim(1)-1]-ct.resolution.y/2;
    dij.ctGrid.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
else
    dij.ctGrid.x = ct.x;
    dij.ctGrid.y = ct.y;   
    dij.ctGrid.z = ct.z;
end

dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);


%Create Dose Grid
dij.doseGrid = this.doseGrid;

%One can provide the dose grid directly (in the future, variable grids would be possible with this)
if ~all(isfield(dij.doseGrid,{'x','y','z'}))
    dij.doseGrid.x = dij.ctGrid.x(1):this.doseGrid.resolution.x:dij.ctGrid.x(end);
    dij.doseGrid.y = dij.ctGrid.y(1):this.doseGrid.resolution.y:dij.ctGrid.y(end);
    dij.doseGrid.z = dij.ctGrid.z(1):this.doseGrid.resolution.z:dij.ctGrid.z(end);
end

dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);
matRad_cfg.dispInfo('Dose grid has dimensions %dx%dx%d\n',dij.doseGrid.dimensions(1),dij.doseGrid.dimensions(2),dij.doseGrid.dimensions(3));

dij.doseGrid.isoCenterOffset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
    dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
    dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];

%TODO: Maybe we should not do this in the preprocessing if it allows us to not
%change the stf
for i = 1:numel(stf)    
    stf(i).isoCenter = stf(i).isoCenter + dij.doseGrid.isoCenterOffset;
end



% meta information for dij
dij.numOfBeams         = numel(stf);
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if this.calcDoseDirect
    this.numOfColumnsDij      = length(stf);
else
    this.numOfColumnsDij      = dij.totalNumOfBixels;
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(this.numOfColumnsDij,1);
dij.rayNum   = NaN*ones(this.numOfColumnsDij,1);
dij.beamNum  = NaN*ones(this.numOfColumnsDij,1);

%Default MU calibration
dij.minMU               = zeros(this.numOfColumnsDij,1);
dij.maxMU               = inf(this.numOfColumnsDij,1);
dij.numOfParticlesPerMU = 1e6*ones(this.numOfColumnsDij,1);


% Allocate space for dij.physicalDose sparse matrix, assume 1%nnz
for i = 1:dij.numOfScenarios   
    nnzEstimate = floor(0.01*(dij.doseGrid.numOfVoxels*this.numOfColumnsDij));
    dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,nnzEstimate);
end

% take only voxels inside patient
VctGrid = [cst{:,4}];
VctGrid = unique(vertcat(VctGrid{:}));

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
this.VdoseGrid = find(matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));

% save vct grid as own property in order to allow sub-classes
% to access it
this.VctGrid = VctGrid;

% Convert CT subscripts to linear indices.
[this.yCoordsV_vox, this.xCoordsV_vox, this.zCoordsV_vox] = ind2sub(ct.cubeDim,this.VctGrid);


% Convert CT subscripts to coarse linear indices.
[this.yCoordsV_voxDoseGrid, this.xCoordsV_voxDoseGrid, this.zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,this.VdoseGrid);

% load machine file from base data folder
this.machine = this.loadMachine(stf(1).radiationMode,stf(1).machine);

this.doseGrid = dij.doseGrid;

end

