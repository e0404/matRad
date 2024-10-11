function [dij] = initDoseCalc(this,ct,cst,stf)
% matRad_DoseEngine.initDoseCalc: Interface for dose calculation
%   method for setting and preparing the inition parameters for the
%   dose calculation.
%   Should be called at the beginning of calcDose method.
%   Can be expanded or changed by overwriting this method and calling
%   the superclass method inside of it
%
% call:
%   [dij] = matRad_DoseEngine.initDoseCalc(this,ct,cst,stf)
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
    this.hWaitbar = waitbar(0,msg,'Color',matRad_cfg.gui.backgroundColor,'DefaultTextColor',matRad_cfg.gui.textColor);
    matRad_applyThemeToWaitbar(this.hWaitbar);
    % prevent closure of waitbar and show busy state
    set(this.hWaitbar,'pointer','watch');
end
this.lastProgressUpdate = tic;

machine = unique({stf.machine});
radiationMode = unique({stf.radiationMode});
if numel(machine) ~= 1 || numel(radiationMode) ~= 1
    matRad_cfg.dispError('machine and radiation mode need to be unique within supplied stf!');
end
%extract strings from cell
machine = machine{1};
radiationMode = radiationMode{1};

%Scenario Model
if ~isa(this.multScen,'matRad_ScenarioModel')
    this.multScen = matRad_multScen(ct,this.multScen);
end

% load machine file from base data folder
this.machine = this.loadMachine(radiationMode,machine);

%Biological Model
if ~isa(this.bioModel,'matRad_BiologicalModel')
    this.bioModel = matRad_BiologicalModel.validate(this.bioModel,radiationMode, this.providedQuantities(this.machine));
end

if any(strcmp(this.bioModel.requiredQuantities, 'LET'))

    this.calcLET = true;
end

dij = struct();

if matRad_ispropCompat(this.bioModel, 'RBE') && ~isnan(this.bioModel.RBE)
    dij.RBE = this.bioModel.RBE; 
end

%store CT grid
dij.ctGrid.resolution = ct.resolution;

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
ct = matRad_getWorldAxes(ct);

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;   
dij.ctGrid.z = ct.z;

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

dij.doseGrid.cubeCoordOffset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
    dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
    dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];

% meta information for dij
dij.numOfBeams         = numel(stf);
dij.numOfScenarios     = this.multScen.totNumScen;

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

if isempty(this.voxelSubIx)
    % take only voxels inside patient
    tmpVctGridScen = cell(1,ct.numOfCtScen);
    for s = 1:ct.numOfCtScen
        tmpScen = cellfun(@(c) c{s},cst(:,4),'UniformOutput',false);
        tmpVctGridScen{s} = unique(vertcat(tmpScen{:}));    
    end
else
    if iscell(this.voxelSubIx)
        tmpVctGridScen = cell(1,ct.numOfCtScen);
        for s = 1:ct.numOfCtScen
            tmpVctGridScen{s} = this.voxelSubIx;    
        end
    else
        tmpVctGridScen = this.voxelSubIx;
    end    
end
this.VctGrid = unique(vertcat(tmpVctGridScen{:}));
% No we find the subindexes for the indivdual scenarios. This helps us
% doing a subselection later on.
this.VctGridScenIx = cellfun(@(c) ismember(this.VctGrid,c),tmpVctGridScen,'UniformOutput',false);


tmpVdoseGridScen = cell(1,ct.numOfCtScen);
for s = 1:ct.numOfCtScen
    % receive linear indices and grid locations from the dose grid
    tmpCube    = zeros(ct.cubeDim);
    tmpCube(tmpVctGridScen{s}) = 1;
    
    % interpolate cube
    tmpVdoseGridScen{s} = find(matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
        dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));
end
this.VdoseGrid = unique(vertcat(tmpVdoseGridScen{:}));
this.VdoseGridScenIx = cellfun(@(c) ismember(this.VdoseGrid,c), tmpVdoseGridScen,'UniformOutput',false);


% Convert CT subscripts to world coordinates.
this.voxWorldCoords = matRad_cubeIndex2worldCoords(this.VctGrid,dij.ctGrid);

% Convert dosegrid subscripts to world coordinates
this.voxWorldCoordsDoseGrid = matRad_cubeIndex2worldCoords(this.VdoseGrid,dij.doseGrid);

%Create helper masks
this.VdoseGridMask = false(dij.doseGrid.numOfVoxels,1);
this.VdoseGridMask(this.VdoseGrid) = true;

this.VctGridMask = false(prod(ct.cubeDim),1);
this.VctGridMask(this.VctGrid) = true;

this.doseGrid = dij.doseGrid;

%Voxel selection for dose calculation
% ser overlap prioriites
this.cstDoseGrid = matRad_setOverlapPriorities(cst);

% resizing cst to dose cube resolution
this.cstDoseGrid = matRad_resizeCstToGrid(this.cstDoseGrid,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
   dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

%structures that are selected here will be included in dose calculation over the robust scenarios
this.robustVoxelsOnGrid = matRad_selectVoxelsFromCst(this.cstDoseGrid, dij.doseGrid, this.selectVoxelsInScenarios);

end

