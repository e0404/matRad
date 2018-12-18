function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ompMC monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseMc(ct,stf,pln,cst,visBool)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   visBool:                    binary switch to enable visualization
% output
%   dij:                        matRad dij struct
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disable visualiazation by default
if nargin < 5
    visBool = false;
end

%%
if exist('matRad_ompInterface') ~= 3
    try
        eval(['mex -largeArrayDims ' fileparts(mfilename('fullpath')) filesep 'submodules' filesep 'ompMC' filesep 'matRad_ompInterface.c']);
    catch
        error('Could not find/generate mex interface for MC dose calculation');
    end
end


% to guarantee downwards compatibility with data that does not have
% ct.x/y/z
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[1:ct.cubeDim(1)];
    ct.y = ct.resolution.y*[1:ct.cubeDim(2)];
    ct.z = ct.resolution.z*[1:ct.cubeDim(3)];
end

% set grids
if ~isfield(pln,'propDoseCalc') || ...
   ~isfield(pln.propDoseCalc,'doseGrid') || ...
   ~isfield(pln.propDoseCalc.doseGrid,'resolution')
    % default values
    dij.doseGrid.resolution.x = 2.5; % [mm]
    dij.doseGrid.resolution.y = 2.5; % [mm]
    dij.doseGrid.resolution.z = 2.5;   % [mm]
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

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

% take only voxels inside patient
VctGrid = [cst{:,4}];
VctGrid = unique(vertcat(VctGrid{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(dij.ctGrid.numOfVoxels,1);
eraseCtDensMask(VctGrid) = 0;
for i = 1:ct.numOfCtScen
    ct.cubeHU{i}(eraseCtDensMask == 1) = -1024;
end

% downsample ct
for s = 1:dij.numOfScenarios
    HUcube{s} =  interp3(dij.ctGrid.y,  dij.ctGrid.x',  dij.ctGrid.z,ct.cubeHU{s}, ...
                         dij.doseGrid.y,dij.doseGrid.x',dij.doseGrid.z,'linear');
end

%% Setup OmpMC options / parameters
ompMCgeo.ctRes = [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z];

%display options
ompMCoptions.verbose = true;

% start MC control          
ompMCoptions.nHistories = 5000;
ompMCoptions.nBatches = 10;
ompMCoptions.randomSeeds = [97 33];

%start source definition      
ompMCoptions.spectrumFile = [pwd filesep 'submodules' filesep 'ompMC' filesep 'spectra' filesep 'mohan6.spectrum'];
ompMCoptions.monoEnergy = 0.1; 
ompMCoptions.charge = 0;
ompMCoptions.colliBounds = [22.5 27.5 22.5 27.5];
ompMCoptions.ssd = 90.0; %This has to be calculated by matRad?
                                                                    
% start MC transport
ompMCoptions.dataFolder = [pwd filesep 'submodules' filesep 'ompMC' filesep 'data' filesep];
ompMCoptions.pegsFile = [pwd filesep 'submodules' filesep 'ompMC' filesep 'pegs4' filesep '700icru.pegs4dat'];
ompMCoptions.pgs4formFile = [pwd filesep 'submodules' filesep 'ompMC' filesep 'pegs4' filesep 'pgs4form.dat'];

ompMCoptions.global_ecut = 0.700;
ompMCoptions.global_pcut = 0.010; 

% Relative Threshold for dose
ompMCoptions.relDoseThreshold = 0.01;

% Output folders
ompMCoptions.outputFolder = [pwd filesep 'submodules' filesep 'ompMC' filesep 'output'];

% Create Material Density Cube
materialFile = [pwd filesep 'submodules' filesep 'ompMC' filesep 'data' filesep '700icru.pegs4dat'];
material = cell(4,5);
material{1,1} = 'AIR700ICRU';
material{1,2} = -1024; 
material{1,3} = -974;
material{1,4} = 0.001;
material{1,5} = 0.044;
material{2,1} = 'LUNG700ICRU';
material{2,2} = -974; 
material{2,3} = -724;
material{2,4} = 0.044; 
material{2,5} = 0.302;
material{3,1} = 'ICRUTISSUE700ICRU';
material{3,2} = -724; 
material{3,3} = 101;
material{3,4} = 0.302; 
material{3,5} = 1.101;
material{4,1} = 'ICRPBONE700ICRU';
material{4,2} = 101; 
material{4,3} = 1976;
material{4,4} = 1.101; 
material{4,5} = 2.088;

% conversion from HU to densities & materials
for s = 1:dij.numOfScenarios

    % projecting out of bounds HU values where necessary
    if max(HUcube{s}(:)) > material{end,3}
        warning('projecting out of range HU values');
        HUcube{s}(HUcube{s}(:) > material{end,3}) = material{end,3};
    end
    if min(HUcube{s}(:)) < material{1,2}
        warning('projecting out of range HU values');
        HUcube{s}(HUcube{s}(:) < material{1,2}) = material{1,2};
    end

    % find material index
    cubeMatIx{s} = NaN*ones(dij.doseGrid.dimensions,'int32');
    for i = size(material,1):-1:1
        cubeMatIx{s}(HUcube{s} <= material{i,3}) = i;
    end
    
    % create an artificial HU lookup table
    hlut = [];
    for i = 1:size(material,1)       
        hlut = [hlut;material{i,2} material{i,4};material{i,3}-1e-10 material{i,5}]; % add eps for interpolation
    end
    
    cubeRho{s} = interp1(hlut(:,1),hlut(:,2),HUcube{s});

end

ompMCgeo.material = material;
ompMCgeo.materialFile = materialFile;

scale = 10; % to convert to cm

ompMCgeo.xBounds = (dij.doseGrid.resolution.x * (0.5 + [0:dij.doseGrid.dimensions(1)])) ./ scale;
ompMCgeo.yBounds = (dij.doseGrid.resolution.y * (0.5 + [0:dij.doseGrid.dimensions(2)])) ./ scale;
ompMCgeo.zBounds = (dij.doseGrid.resolution.z * (0.5 + [0:dij.doseGrid.dimensions(3)])) ./ scale;

%% debug visualization
if visBool
    
    figure
    hold on

    axis equal
    
    % ct box
    ctCorner1 = [ompMCgeo.xBounds(1) ompMCgeo.yBounds(1) ompMCgeo.zBounds(1)];
    ctCorner2 = [ompMCgeo.xBounds(end) ompMCgeo.yBounds(end) ompMCgeo.zBounds(end)];
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner1(3) ctCorner1(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner2(2)],[ctCorner2(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner1(2) ctCorner1(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner1(1) ctCorner1(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )
    plot3([ctCorner2(1) ctCorner2(1)],[ctCorner2(2) ctCorner2(2)],[ctCorner1(3) ctCorner2(3)],'k' )
        
    xlabel('x [cm]')
    ylabel('y [cm]')
    zlabel('z [cm]')

    rotate3d on
    
end

%% Create beamlet source
numOfBixels = zeros(dij.numOfBeams, 1);
beamSource = zeros(dij.numOfBeams, 3);

bixelCorner = zeros(dij.totalNumOfBixels,3);
bixelSide1 = zeros(dij.totalNumOfBixels,3);
bixelSide2 = zeros(dij.totalNumOfBixels,3);

counter = 0;

for i = 1:dij.numOfBeams % loop over all beams
   
    % define beam source in physical coordinate system in cm
    beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;

    numOfBixels(i) = stf(i).numOfRays;

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
        % get bixel corner and delimiting vectors.
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        currCorner = (stf(i).ray(j).beamletCornersAtIso(1,:) + stf(i).isoCenter) ./ scale;
        bixelCorner((i-1)*stf(i).numOfRays + j,:) = currCorner;
        bixelSide1((i-1)*stf(i).numOfRays + j,:) = (stf(i).ray(j).beamletCornersAtIso(2,:) + stf(i).isoCenter) ./ scale - currCorner;
        bixelSide2((i-1)*stf(i).numOfRays + j,:) = (stf(i).ray(j).beamletCornersAtIso(4,:) + stf(i).isoCenter) ./ scale - currCorner;
        
        if visBool
            for k = 1:4
                currCornerVis = (stf(i).ray(j).beamletCornersAtIso(k,:) + stf(i).isoCenter)/10;
                % rays connecting source and ray corner
                plot3([beamSource(i,1) currCornerVis(1)],[beamSource(i,2) currCornerVis(2)],[beamSource(i,3) currCornerVis(3)],'y')
                % connection between corners
                lRayCorner = (stf(i).ray(j).beamletCornersAtIso(mod(k,4) + 1,:) + stf(i).isoCenter)/10;
                plot3([lRayCorner(1) currCornerVis(1)],[lRayCorner(2) currCornerVis(2)],[lRayCorner(3) currCornerVis(3)],'r')
            end
        end
        
    end
        
end

ompMCsource.nBeams = dij.numOfBeams;
ompMCsource.xSource = beamSource(:,1);
ompMCsource.ySource = beamSource(:,2);
ompMCsource.zSource = beamSource(:,3);

ompMCsource.nBixels = numOfBixels(:);
ompMCsource.xCorner = bixelCorner(:,1);
ompMCsource.yCorner = bixelCorner(:,2);
ompMCsource.zCorner = bixelCorner(:,3);

ompMCsource.xSide1 = bixelSide1(:,1);
ompMCsource.ySide1 = bixelSide1(:,2);
ompMCsource.zSide1 = bixelSide1(:,3);

ompMCsource.xSide2 = bixelSide2(:,1);
ompMCsource.ySide2 = bixelSide2(:,2);
ompMCsource.zSide2 = bixelSide2(:,3);

if visBool
    plot3(ompMCsource.xSource,ompMCsource.zSource,ompMCsource.ySource,'rx')
end

%% Call the OmpMC interface
%initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for photons...');
% show busy state
set(figureWait,'pointer','watch');

fprintf('matRad: OmpMC photon dose calculation... ');

%run over all scenarios
for s = 1:dij.numOfScenarios
    ompMCgeo.isoCenter = [stf(:).isoCenter];
    dij.physicalDose{s} = matRad_ompInterface(cubeRho{s},cubeMatIx{s},ompMCgeo,ompMCsource,ompMCoptions);
end

fprintf('matRad: done!\n');

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

end
