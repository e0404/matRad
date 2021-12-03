matRad_cfg =  MatRad_Config.instance();

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[0:ct.cubeDim(2)-1]-ct.resolution.x/2;
    ct.y = ct.resolution.y*[0:ct.cubeDim(1)-1]-ct.resolution.y/2;
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
   ~isfield(pln.propDoseCalc,'doseGrid') || ...
   ~isfield(pln.propDoseCalc.doseGrid,'resolution')
    % default values
    dij.doseGrid.resolution = matRad_cfg.propDoseCalc.defaultResolution;
else
    % take values from pln strcut
    dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
    dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
    dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
end

dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

dij.doseGrid.dimensions  = [numel(dij.doseGrid.y) numel(dij.doseGrid.x) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.y) numel(dij.ctGrid.x) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

% adjust isocenter internally for different dose grid
offset = [dij.doseGrid.resolution.x - dij.ctGrid.resolution.x ...
          dij.doseGrid.resolution.y - dij.ctGrid.resolution.y ...
          dij.doseGrid.resolution.z - dij.ctGrid.resolution.z];
   
% meta information for dij
switch pln.radiationMode
    case {'photons','protons','carbon'}

        dij.numOfBeams          = pln.propStf.numOfBeams;
        dij.numOfScenarios      = 1;
        dij.numOfRaysPerBeam    = [stf(:).numOfRays];
        dij.totalNumOfBixels    = sum([stf(:).totalNumOfBixels]);
        dij.totalNumOfRays      = sum(dij.numOfRaysPerBeam);
    case 'brachy'
        dij.numOfScenarios      = 1;
        dij.numOfBeams          = 1;
        dij.beamNum             = 1;
        dij.numOfNeedles        = stf.numOfNeedles;
        dij.numOfSeedsPerNeedle = stf.numOfSeedsPerNeedle;
        dij.totalNumOfSeeds     = dij.numOfNeedles*dij.numOfSeedsPerNeedle;
        dij.totalNumOfBixels    = dij.totalNumOfSeeds;
        
end


% Allocate space for dij.physicalDose sparse matrix
dij.physicalDose = spalloc(dij.doseGrid.numOfVoxels,dij.totalNumOfBixels,1);


% take only voxels inside patient
VctGrid = [cst{:,4}];
VctGrid = unique(vertcat(VctGrid{:}));

% ignore densities outside of contours
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'ignoreOutsideDensities')
    ignoreOutsideDensities = matRad_cfg.propDoseCalc.defaultIgnoreOutsideDensities;
else
    ignoreOutsideDensities = pln.propDoseCalc.ignoreOutsideDensities;
end

if ignoreOutsideDensities
    eraseCtDensMask = ones(prod(ct.cubeDim),1);
    eraseCtDensMask(VctGrid) = 0;
    for i = 1:ct.numOfCtScen
        ct.cube{i}(eraseCtDensMask == 1) = 0;
    end
end

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,VctGrid);

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
VdoseGrid = find(matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z,tmpCube, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest'));

% Convert CT subscripts to coarse linear indices.
[yCoordsV_voxDoseGrid, xCoordsV_voxDoseGrid, zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,VdoseGrid);

% load base data% load machine file
fileName = [pln.radiationMode '_' pln.machine '.mat'];
try
   load(fullfile(matRad_cfg.matRadRoot,'basedata',fileName));
catch
   matRad_cfg.dispError('Could not find the following machine file: %s\n',fileName); 
end
dij.basedata = machine;
% compute SSDs
% stf = matRad_computeSSD(stf,ct); 
