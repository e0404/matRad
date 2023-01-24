function [dij,ct,cst,stf,pln] = calcDoseInit(this,ct,cst,stf)
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

%   stf:            matRad stf struct
%   pln:            matRad pln struct
%
% returns:
%   ct:             matRad ct  struct
%   stf:            matRad stf struct
%   pln:            matRad pln struct
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

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1]-ct.resolution.x/2;
    ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1]-ct.resolution.y/2;
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
end

dij = struct();

%One can provide the dose grid directly (in the future, variable grids would be possible with this)
if ~all(isfield(this.doseGrid,{'x','y','z'}))
    dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
    dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
    dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);
end

dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

%store CT grid
dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

dij.doseGrid.isoCenterOffset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
    dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
    dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];

%Maybe we should not do this in the preprocessing if it allows us to not
%change the stf
for i = 1:numel(stf)
    stf(i).isoCenter = stf(i).isoCenter + dij.doseGrid.isoCenterOffset;
end

%If we want to omit HU conversion check if we have a ct.cube ready
if this.useGivenEqDensityCube && ~isfield(ct,'cube')
    matRad_cfg.dispWarning('HU Conversion requested to be omitted but no ct.cube exists! Will override and do the conversion anyway!');
    this.useGivenEqDensityCube = false;
end

% calculate rED or rSP from HU
% Maybe we can avoid duplicating the CT here?
if this.useGivenEqDensityCube
    matRad_cfg.dispInfo('Omitting HU to rED/rSP conversion and using existing ct.cube!\n');
else
    ct = matRad_calcWaterEqD(ct, pln);
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
if this.calcDoseDirect
    this.numOfColumnsDij      = length(stf);
    this.numOfBixelsContainer = 1;
else
    this.numOfColumnsDij      = dij.totalNumOfBixels;
    this.numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
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
    nnzEstimate = 0.01*(dij.doseGrid.numOfVoxels*this.numOfColumnsDij);
    dij.physicalDose{i} = spalloc(dij.doseGrid.numOfVoxels,this.numOfColumnsDij,nnzEstimate);
end

% take only voxels inside patient
VctGrid = [cst{:,4}];
VctGrid = unique(vertcat(VctGrid{:}));

% ignore densities outside of contours
if this.ignoreOutsideDensities
    eraseCtDensMask = ones(prod(ct.cubeDim),1);
    eraseCtDensMask(VctGrid) = 0;
    for i = 1:ct.numOfCtScen
        ct.cube{i}(eraseCtDensMask == 1) = 0;
    end
end

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
this.machine = matRad_loadMachine(pln);

% compute SSDs
stf = matRad_computeSSD(stf,ct,'densityThreshold',this.ssdDensityThreshold);

end

