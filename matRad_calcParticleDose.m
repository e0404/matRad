function dij = matRad_calcParticleDose(ct,stf,pln,cst,param,heteroCorrBio)
% matRad particle dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   param:          (optional) structure defining additional parameter
%                   param.calcDoseDirect boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('param','var')
    if ~isfield(param,'logLevel')
        param.logLevel = 1;
    end
    % default: dose influence matrix computation
    if ~isfield(param,'calcDoseDirect')
        param.calcDoseDirect = false;
    end
else
    param.calcDoseDirect = false;
    param.subIx          = [];
    param.logLevel       = 1;
end

if ~exist('heteroCorrBio','var') || isempty(heteroCorrBio)
   heteroCorrBio = false;
end

% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[0:ct.cubeDim(1)-1]-ct.resolution.x/2;
    ct.y = ct.resolution.y*[0:ct.cubeDim(2)-1]-ct.resolution.y/2;
    ct.z = ct.resolution.z*[0:ct.cubeDim(3)-1]-ct.resolution.z/2;
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
        ~isfield(pln.propDoseCalc,'doseGrid') || ...
        ~isfield(pln.propDoseCalc.doseGrid,'resolution')
    % default values
    dij.doseGrid.resolution.x = 2.5; % [mm]
    dij.doseGrid.resolution.y = 2.5; % [mm]
    dij.doseGrid.resolution.z = 3;   % [mm]
else
    % take values from pln strcut
    dij.doseGrid.resolution.x = pln.propDoseCalc.doseGrid.resolution.x;
    dij.doseGrid.resolution.y = pln.propDoseCalc.doseGrid.resolution.y;
    dij.doseGrid.resolution.z = pln.propDoseCalc.doseGrid.resolution.z;
end

dij.doseGrid.x = ct.x(1):dij.doseGrid.resolution.x:ct.x(end);
dij.doseGrid.y = ct.y(1):dij.doseGrid.resolution.y:ct.y(end);
dij.doseGrid.z = ct.z(1):dij.doseGrid.resolution.z:ct.z(end);

dij.doseGrid.dimensions  = [numel(dij.doseGrid.x) numel(dij.doseGrid.y) numel(dij.doseGrid.z)];
dij.doseGrid.numOfVoxels = prod(dij.doseGrid.dimensions);

dij.ctGrid.resolution.x = ct.resolution.x;
dij.ctGrid.resolution.y = ct.resolution.y;
dij.ctGrid.resolution.z = ct.resolution.z;

dij.ctGrid.x = ct.x;
dij.ctGrid.y = ct.y;
dij.ctGrid.z = ct.z;

dij.ctGrid.dimensions  = [numel(dij.ctGrid.x) numel(dij.ctGrid.y) numel(dij.ctGrid.z)];
dij.ctGrid.numOfVoxels = prod(dij.ctGrid.dimensions);

if param.logLevel == 1
    % initialize waitbar
    figureWait = waitbar(0,'calculate dose influence matrix for particles...');
    % prevent closure of waitbar and show busy state
    set(figureWait,'pointer','watch');
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln, param);

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = pln.multScen.numOfCtScen;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

if param.calcDoseDirect
    numOfColumnsDij      = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij      = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

% Allocate space for dij.physicalDose sparse matrix
for ctScen = 1:pln.multScen.numOfCtScen
    for shiftScen = 1:pln.multScen.totNumShiftScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                dij.physicalDose{ctScen,shiftScen,rangeShiftScen} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
            end
            
        end
    end
end


% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
doseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);


% if biological optimization considering a variable RBE is true then create alphaDose and betaDose containers and sparse matrices
if pln.bioParam.bioOpt
    
    alphaDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for shiftScen = 1:pln.multScen.totNumShiftScen
            for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                
                if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                    dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}        = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                    dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}     = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                end
                
            end
            
        end
        
    end
    
end


% Only take voxels inside patient at full ct resolution
if isfield(param,'subIx') && ~isempty(param.subIx)
    VctGrid = param.subIx;
else
    VctGrid = [cst{:,4}];
    VctGrid = unique(vertcat(VctGrid{:}));
end

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(VctGrid) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,VctGrid);

% receive linear indices and grid locations from the dose grid
tmpCube    = zeros(ct.cubeDim);
tmpCube(VctGrid) = 1;
% interpolate cube
VdoseGrid = find(interp3(dij.ctGrid.y,dij.ctGrid.x,   dij.ctGrid.z,tmpCube, ...
    dij.doseGrid.y,dij.doseGrid.x',dij.doseGrid.z,'nearest'));

% Convert CT subscripts to coarse linear indices.
[yCoordsV_voxDoseGrid, xCoordsV_voxDoseGrid, zCoordsV_voxDoseGrid] = ind2sub(dij.doseGrid.dimensions,VdoseGrid);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
    load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
    matRad_dispToConsole(['Could not find the following machine file: ' fileName ],param,'error');
end

% check if we do heterogeneity correction
matRad_dispToConsole('check if we do heterogeneity correction. \n',param,'info');
calcHeteroCorr = false;
for j = 1:size(cst,1)
    if isfield(cst{j,5},'HeterogeneityCorrection')
        if strcmp(cst{j,5}.HeterogeneityCorrection,'Lung')
            calcHeteroCorr = true;
            break;
        end
    end
end

if calcHeteroCorr
    calcHeteroCorrStruct.cube = {zeros(ct.cubeDim)};
    calcHeteroCorrStruct.cubeDim = ct.cubeDim;
    calcHeteroCorrStruct.numOfCtScen = pln.multScen.numOfCtScen;
    calcHeteroCorrStruct.resolution = ct.resolution;
    
    % book keeping - this is necessary since pln is not used in optimization or
    % matRad_calcCubes
    
    for j = 1:size(cst,1)
        if isfield(cst{j,5},'HeterogeneityCorrection')
            if strcmp(cst{j,5}.HeterogeneityCorrection,'Lung')
                calcHeteroCorrStruct.cube{1}(cst{j,4}{1}) = ct.cube{1}(cst{j,4}{1});
            else
                error(['No heterogeneity correction method implemented for ' ...
                    cst{j,5}.HeterogeneityCorrection '.']);
            end
        end
    end
end

if isfield(pln,'propDoseCalc') && ...
        isfield(pln.propDoseCalc,'calcLET') && ...
        pln.propDoseCalc.calcLET
    if isfield(machine.data,'LET')
        
        letDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
        
        for ctScen = 1:pln.multScen.numOfCtScen
            for shiftScen = 1:pln.multScen.totNumShiftScen
                for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                    
                    if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                        dij.mLETDose{ctScen,shiftScen,rangeShiftScen} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                    end
                    
                end
            end
        end
        
    else
        matRad_dispToConsole('LET not available in the machine data. LET will not be calculated.',param,'warning');
    end
end

if strcmp(pln.bioParam.model,'constRBE')
    dij.RBE = pln.bioParam.RBE;
end

% ser overlap prioriites
cst = matRad_setOverlapPriorities(cst);

% resizing cst to dose cube resolution
cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
    dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);

vTissueIndex = zeros(size(VdoseGrid,1),1);
% generates tissue class matrix for biological treatment planning and alpha_x, beta_x, vectors
if pln.bioParam.bioOpt
    
    dij.alphaX  = zeros(dij.doseGrid.numOfVoxels,1);
    dij.betaX   = zeros(dij.doseGrid.numOfVoxels,1);
    dij.abX     = zeros(dij.doseGrid.numOfVoxels,1);  % alpha beta ratio
    
    matRad_dispToConsole(['matRad: Please use same priorities and tissue classes for optimization \n'],param,'info');
    %set overlap priorities
    cst = matRad_setOverlapPriorities(cst);
    
    % create radiosensitivity vectors
    for i = 1:size(cst,1)
        % check if cst is compatible
        if ~isfield(cst{i,5},'alphaX') || ~isfield(cst{i,5},'betaX')
            cst{i,5}.alphaX = 0.1;
            cst{i,5}.betaX  = 0.05;
            matRad_dispToConsole(['matRad: using default alpha_x and beta_x parameters for ' cst{i,2} ' \n'],param,'warning');
        end
        
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
            dij.alphaX(cst{i,4}{1}) = cst{i,5}.alphaX;
            dij.betaX(cst{i,4}{1})  = cst{i,5}.betaX;
        end
        
    end
    
    
    
    
    % retrieve photon LQM parameter for the current dose grid voxels
    [dij.alphaX,dij.betaX] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);
    
    
    
    dij.abX(dij.betaX>0) = dij.alphaX(dij.betaX>0)./dij.betaX(dij.betaX>0);
    
    % only if LEM is used corresponding bio data must be available in the base data set
    if strcmp(pln.radiationMode,'carbon') || strcmp(pln.radiationMode,'helium')
        
        for i = 1:size(cst,1)
            [~, row] = ismember(vertcat(cst{i,4}{:}),VdoseGrid,'rows');
            % check if cst is compatiable
            if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
                
                if ~isempty(cst{i,5})
                    
                    IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                        ismember(machine.data(1).betaX, cst{i,5}.betaX));
                    
                    % check consistency of biological baseData and cst settings
                    if ~isempty(IdxTissue)
                        vTissueIndex(row) = IdxTissue;
                    else
                        matRad_dispToConsole('biological base data and cst inconsistent \n',param,'error');
                    end
                    
                else
                    vTissueIndex(row) = 1;
                    matRad_dispToConsole(['matRad: tissue type of ' cst{i,2} ' was set to 1  \n'],param,'info');
                end
            else
                matRad_dispToConsole('base data is incomplement - alphaX and/or betaX is missing',param,'error');
            end
        end
        
        matRad_dispToConsole('done. \n',param,'info');
        
    end
    
end

dij.vTissueIndex = vTissueIndex;

ctScen  = 1;        % current ct scenario
matRad_dispToConsole('matRad: Particle dose calculation... \n',param,'info');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over all shift scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for shiftScen = 1:pln.multScen.totNumShiftScen
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end
    
    matRad_dispToConsole(['shift scenario ' num2str(shiftScen) ' of ' num2str(pln.multScen.totNumShiftScen) ':  \n'],param,'info');
    matRad_dispToConsole('matRad: Particle dose calculation... \n',param,'info');
    
    counter = 0;
    
    % compute SSDs only for first scenario
    stf = matRad_computeSSD(stf,ct,param);
    
    for i = 1:numel(stf) % loop over all beams
        
        matRad_dispToConsole(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ':  \n'],param,'info');
        
        % remember beam and bixel number
        if param.calcDoseDirect
            dij.beamNum(i)    = i;
            dij.rayNum(i)     = i;
            dij.bixelNum(i)   = i;
        end
        
        bixelsPerBeam = 0;
        
        % convert voxel indices to real coordinates using iso center of beam i
        xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
        yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
        zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
        coordsV  = [xCoordsV yCoordsV zCoordsV];
        
        xCoordsVdoseGrid = xCoordsV_voxDoseGrid(:) * dij.doseGrid.resolution.x-stf(i).isoCenter(1);
        yCoordsVdoseGrid = yCoordsV_voxDoseGrid(:) * dij.doseGrid.resolution.y-stf(i).isoCenter(2);
        zCoordsVdoseGrid = zCoordsV_voxDoseGrid(:) * dij.doseGrid.resolution.z-stf(i).isoCenter(3);
        coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];
        
        % Get Rotation Matrix
        % Do not transpose matrix since we usage of row vectors &
        % transformation of the coordinate system need double transpose
        
        % rotation around Z axis (gantry)
        rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);
        
        % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
        rot_coordsV         = coordsV*rotMat_system_T;
        rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;
        
        rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
        rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
        rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);
        
        rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
        rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
        rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);
        
        
        % Calcualte radiological depth cube
        lateralCutoffRayTracing = 50;
        matRad_dispToConsole('matRad: calculate radiological depth cube...',param,'info');
        radDepthVctGrid = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,lateralCutoffRayTracing);
        matRad_dispToConsole('done. \n',param,'info');
        
        if calcHeteroCorr
            fprintf('matRad: calculate radiological depth cube for heterogeneity correction...');
            heteroCorrDepthV = matRad_rayTracing(stf(i),calcHeteroCorrStruct,VctGrid,rot_coordsV,lateralCutoffRayTracing);
            % HETERO interpolate hetero depth cube to dose grid resolution
            heteroCorrDepthV = matRad_interpRadDepth...
                (ct,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,heteroCorrDepthV);
            fprintf('done.\n');
        end
        
        % interpolate radiological depth cube to dose grid resolution
        radDepthVdoseGrid = matRad_interpRadDepth...
            (ct,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
        
        
        
        % limit rotated coordinates to positions where ray tracing is availabe
        rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);
        
        % Determine lateral cutoff
        matRad_dispToConsole('matRad: calculate lateral cutoff...',param,'info');
        cutOffLevel          = 0.99;
        visBoolLateralCutOff = 0;
        % compute cut off only for the first scenario
        machine              = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
        matRad_dispToConsole('done. \n',param,'info');
        
        for j = 1:stf(i).numOfRays % loop over all rays
            
            if ~isempty(stf(i).ray(j).energy)
                
                % find index of maximum used energy (round to keV for numerical reasons
                energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);
                
                maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);
                
                % Ray tracing for beam i and ray j
                [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                    stf(i).sourcePoint_bev, ...
                    stf(i).ray(j).targetPoint_bev, ...
                    machine.meta.SAD, ...
                    find(~isnan(radDepthVdoseGrid{1})), ...
                    maxLateralCutoffDoseCalc);
                               
                if calcHeteroCorr
                    heteroCorrDepths = heteroCorrDepthV{1}(ix);
                end
                
                % just use tissue classes of voxels found by ray tracer
  %              if pln.bioParam.bioOpt
                    vTissueIndex_j = vTissueIndex(ix,:);
  %              end
                
                for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
                    
                    counter       = counter + 1;
                    bixelsPerBeam = bixelsPerBeam + 1;
                    
                    if param.logLevel <= 2
                        % Display progress and update text only 200 times
                        if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                            matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                        end
                        if param.logLevel == 1
                            % update waitbar only 100 times if it is not closed
                            if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                                waitbar(counter/dij.totalNumOfBixels,figureWait);
                            end
                        end
                    end
                    
                    % remember beam and bixel number
                    if ~param.calcDoseDirect
                        dij.beamNum(counter)  = i;
                        dij.rayNum(counter)   = j;
                        dij.bixelNum(counter) = k;
                    end
                    
                    % find energy index in base data
                    energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                    
                    % create offset vector to account for additional offsets modelled in the base data and a potential
                    % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                    % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
                    % interpolations in matRad_calcParticleDoseBixel() and avoid extrapolations.
                    offsetRadDepth = machine.data(energyIx).offset - stf(i).ray(j).rangeShifter(k).eqThickness;
                    
                    for ctScen = 1:pln.multScen.numOfCtScen
                        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                            
                            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                
                                radDepths = radDepthVdoseGrid{ctScen}(ix);
                                
                                % manipulate radDepthCube for range scenarios
                                if pln.multScen.relRangeShift(rangeShiftScen) ~= 0 || pln.multScen.absRangeShift(rangeShiftScen) ~= 0
                                    radDepths = radDepths +...                                                       % original cube
                                        radDepthVdoseGrid{ctScen}(ix)*pln.multScen.relRangeShift(rangeShiftScen) +... % rel range shift
                                        pln.multScen.absRangeShift(rangeShiftScen);                                   % absolute range shift
                                    radDepths(radDepths < 0) = 0;
                                end
                                
                                % find depth depended lateral cut off
                                if cutOffLevel >= 1
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth;
                                elseif cutOffLevel < 1 && cutOffLevel > 0
                                    % perform rough 2D clipping
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                        radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);
                                    
                                    % peform fine 2D clipping
                                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                            (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
                                    end
                                    
                                    % empty bixels may happen during recalculation of error
                                    % scenarios -> skip to next bixel
                                    if ~any(currIx)
                                        continue;
                                    end
                                    
                                    % adjust radDepth according to range shifter
                                    currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;
                                    if calcHeteroCorr
                                        currHeteroCorrDepths = heteroCorrDepths(currIx);
                                    end
                                    
                                    % calculate initial focus sigma
                                    sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist(stf(i).ray(j).focusIx(k),:)', ...
                                        machine.data(energyIx).initFocus.sigma(stf(i).ray(j).focusIx(k),:)',stf(i).ray(j).SSD);
                                    sigmaIni_sq = sigmaIni^2;
                                    
                                    % consider range shifter for protons if applicable
                                    if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')
                                        
                                        % compute!
                                        sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy, ...
                                            stf(i).ray(j).rangeShifter(k), ...
                                            stf(i).ray(j).SSD{ctScen});
                                        
                                        % add to initial sigma in quadrature
                                        sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;
                                        
                                    end
                                    
                                    % calculate particle dose for bixel k on ray j of beam i
                                    if calcHeteroCorr
                                        bixel = matRad_calcParticleDoseBixel(...
                                            currRadDepths, ...
                                            radialDist_sq(currIx), ...
                                            sigmaIni_sq, ...
                                            machine.data(energyIx), ...
                                            currHeteroCorrDepths, ...
                                            [], ...
                                            heteroCorrBio, vTissueIndex_j(currIx));
                                    else
                                        bixel = matRad_calcParticleDoseBixel(...
                                            currRadDepths, ...
                                            radialDist_sq(currIx), ...
                                            sigmaIni_sq, ...
                                            machine.data(energyIx),[],[],heteroCorrBio, vTissueIndex_j(currIx));
                                    end
                                    
                                    bixelDose = bixel.physDose;
                                    
                                    %bixelDose = bixel.physDose
                                    
                                else
                                    matRad_dispToConsole('cutoff must be a value between 0 and 1',param,'error')
                                end
                                
                                
                                % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                                %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                                %Type               = 'dose';
                                %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);
                                
                                % save dose for every bixel in cell array
                                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);
                                
                                
                                if isfield(dij,'mLETDose')
                                    % calculate particle LET for bixel k on ray j of beam i
                                    depths   = machine.data(energyIx).depths + machine.data(energyIx).offset;
                                    bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,currRadDepths);
                                    bixelLET(isnan(bixelLET)) = 0;
                                    % save LET for every bixel in cell array
                                    letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                end
                                
                                % save alpha_p and beta_p radiosensititvy parameter for every bixel in cell array
                                if pln.bioParam.bioOpt
                                    
                                    if all(isfield(bixel,{'Z_Aij','Z_Bij'}))
                                        bixelAlphaDose =  bixel.L .* bixel.Z_Aij;
                                        bixelBetaDose  =  bixel.L .* bixel.Z_Bij;
                                    else
                                        [bixelAlpha,bixelBeta] = pln.bioParam.calcLQParameter(radDepths(currIx),machine.data(energyIx),vTissueIndex_j(currIx,:),dij.alphaX(VdoseGrid(ix(currIx))),...
                                            dij.betaX(VdoseGrid(ix(currIx))),...
                                            dij.abX(VdoseGrid(ix(currIx))));
                                        bixelAlphaDose =  bixel.physDose .* bixelAlpha;
                                        bixelBetaDose  =  bixel.physDose .* sqrt(bixelBeta);
                                    end
                                    
                                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelAlphaDose,dij.doseGrid.numOfVoxels,1);
                                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}  = sparse(VdoseGrid(ix(currIx)),1,bixelBetaDose,dij.doseGrid.numOfVoxels,1);
                                    
                                end
                                
                                
                            end
                        end
                    end
                    
                    % save computation time and memory by sequentially filling the
                    % sparse matrix dose.dij from the cell array
                    if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                        
                        for ctScen = 1:pln.multScen.numOfCtScen
                            for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                                if ~any(currIx)
                                    continue;
                                end
                                if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                    
                                    if param.calcDoseDirect
                                        if isfield(stf(1).ray(1),'weight') && numel(stf(i).ray(j).weight) >= k
                                            
                                            % score physical dose
                                            dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * doseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                                            
                                            if isfield(dij,'mLETDose') && pln.sampling
                                                % score LETxDose matrices
                                                dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * letDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                                            end
                                            
                                            if pln.bioParam.bioOpt
                                                % score alphaxDose and sqrt(beta)xDose matrices
                                                dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,i)    = dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,i)    + stf(i).ray(j).weight(k) * alphaDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                                                dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * betaDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                                            end
                                        else
                                            matRad_dispToConsole(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)],param,'error');
                                        end
                                    else
                                        
                                        % fill entire dose influence matrix
                                        dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}];
                                        
                                        if isfield(dij,'mLETDose')
                                            % fill entire LETxDose influence matrix
                                            dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [letDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}];
                                        end
                                        
                                        if pln.bioParam.bioOpt
                                            % fill entire alphaxDose influence and sqrt(beta)xDose influence matrices
                                            dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter)    = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}];
                                            dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}];
                                        end
                                        
                                    end  % if direct dose calc
                                end  % shift scen
                            end % range scen
                            
                            
                        end % end % ct scen
                    end
                end % end bixels per ray
                
            end
            
        end %  end ray loop
        
    end % end beam loop
    
    
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
    end
    
end % end shift scenario loop
matRad_dispToConsole('Ended shifting scenario loop. \n',param,'info');
% remove dose influence for voxels outside of segmentations for every ct
% scenario
for i = 1:pln.multScen.numOfCtScen
    
    % generate index set to erase
    tmpIx = [];
    for j = 1:size(cst,1)
        tmpIx = unique([tmpIx; cst{j,4}{i}]);
    end
    ix = setdiff(1:dij.doseGrid.numOfVoxels,tmpIx);
    
    for j = 1:pln.multScen.totNumShiftScen
        for k = 1:pln.multScen.totNumRangeScen
            
            if pln.multScen.scenMask(i,j,k)
                
                dij.physicalDose{i,j,k}(ix,:)      = 0;
                
                if isfield(dij,'mLETDose')
                    dij.mLETDose{i,j,k}(ix,:)      = 0;
                end
                
                if pln.bioParam.bioOpt
                    dij.mAlphaDose{i,j,k}(ix,:)    = 0;
                    dij.mSqrtBetaDose{i,j,k}(ix,:) = 0;
                end
                
            end
            
        end
    end
end


matRad_dispToConsole('calcParticleDose finished. \n',param,'info');

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end
