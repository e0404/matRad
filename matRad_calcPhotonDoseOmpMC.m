function dij = matRad_calcPhotonDoseOmpMC(ct,stf,pln,cst,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad ompMC monte carlo photon dose calculation wrapper
%
% call
%   dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst)
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

% disable visualiazation by default
if nargin < 5
    visBool = false;
end

% meta information for dij
dij.numOfBeams         = pln.propStf.numOfBeams;
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfScenarios     = 1;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

% check if full dose influence data is required
%{
if calcDoseDirect
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end
%}

numOfColumnsDij           = dij.totalNumOfBixels;
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

doseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);

% Allocate space for dij.physicalDose sparse matrix
%for i = 1:dij.numOfScenarios
%    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
%end



% set environment variables for ompMC
%{
if exist(['vmc++' filesep 'bin'],'dir') ~= 7
    error(['Could not locate vmc++ environment. ' ...
          'Please provide the files in the correct folder structure at matRadroot' filesep 'vmc++.']);
else
    VMCPath     = fullfile(pwd , 'vmc++');
    runsPath    = fullfile(VMCPath, 'runs');
    phantomPath = fullfile(VMCPath, 'phantoms');

    setenv('vmc_home',VMCPath);
    setenv('vmc_dir',runsPath);
    setenv('xvmc_dir',VMCPath);
    
    if isunix
        system(['chmod a+x ' VMCPath filesep 'bin' filesep 'vmc_Linux.exe']);
    end
    
end
%}

%% Setup OmpMC options / parameters
%ompMCgeo.nFields = numel(stf);
%ompMCgeo.nBixels = [stf(:).totalNumOfBixels];
ompMCgeo.ctRes = [ct.resolution.x ct.resolution.y ct.resolution.z];

%display options
ompMCoptions.verbose = true;

% start MC control          
ompMCoptions.nHistories = 5000;
ompMCoptions.nBatches = 10;
ompMCoptions.randomSeeds = [97 33];


%start source definition      
ompMCoptions.spectrumFile = [pwd '/ompMC/spectra/mohan6.spectrum'];
ompMCoptions.monoEnergy = 0.1; 
ompMCoptions.charge = 0;
ompMCoptions.colliBounds = [-0.25 0.25 -0.25 0.25];
ompMCoptions.ssd = 90.0; %This has to be calculated by matRad?
                                                                    
% start MC transport
ompMCoptions.dataFolder = [pwd '/ompMC/data/'];
ompMCoptions.pegsFile = [pwd '/ompMC/pegs4/700icru.pegs4dat'];
ompMCoptions.pgs4formFile = [pwd '/ompMC/pegs4/pgs4form.dat'];

ompMCoptions.global_ecut = 0.700;
ompMCoptions.global_pcut = 0.010; 

% Relative Threshold for dose
ompMCoptions.relDoseThreshold = 0.01;

% Output folders
ompMCoptions.outputFolder = [pwd '/ompMC/output/'];


%% Create Material Density Cube
materialFile = [pwd '/ompMC/700icru.pegs4dat'];
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
          
    % find material index
    cubeMatIx{s} = NaN*ones(ct.cubeDim,'int32');
    for i = size(material,1):-1:1
        cubeMatIx{s}(ct.cubeHU{s} <= material{i,3}) = i;
    end
    
    % create an artificial HU lookup table
    hlut = [];
    for i = 1:size(material,1)       
        hlut = [hlut;material{i,2} material{i,4};material{i,3}-1e-10 material{i,5}]; % add eps for interpolation
    end
    
    cubeRho{s} = interp1(hlut(:,1),hlut(:,2),ct.cubeHU{s});

end

ompMCgeo.material = material;
ompMCgeo.materialFile = materialFile;

scale = 10; %cm?

if ~isfield(ct,'x')
    %length = (1:ct.cubeDim(1) * ct.resolution.x) - ct.resolution.x;
    ct.x =  ct.resolution.x * [0:ct.cubeDim(1)-1] - (ct.resolution.x * ct.cubeDim(1))/2;
    %ct.x = (ct.cubeDim(1)/2 - 1:ct.cubeDim(1));
end

if ~isfield(ct,'y')
    ct.y =  ct.resolution.y * [0:ct.cubeDim(2)-1] - (ct.resolution.y * ct.cubeDim(2))/2;
end

if ~isfield(ct,'z')
    ct.z =  ct.resolution.z * [0:ct.cubeDim(3)-1] - (ct.resolution.z * ct.cubeDim(3))/2;
end

ompMCgeo.xBounds = [(ct.x - ct.resolution.x*0.5) (ct.x(ct.cubeDim(1)) + ct.resolution.x)] ./ scale;
ompMCgeo.yBounds = [(ct.y - ct.resolution.y*0.5) (ct.y(ct.cubeDim(2)) + ct.resolution.y)] ./ scale;
ompMCgeo.zBounds = [(ct.z - ct.resolution.z*0.5) (ct.z(ct.cubeDim(3)) + ct.resolution.z)] ./ scale;

%% visualization
if visBool

    clf
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

for i = 1:dij.numOfBeams % loop over all beams
   
    % define beam source in physical coordinate system in cm
    beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;
    %fwrite(fileHandle,beamSource,'double');   
    
    % write number of beamlets into file
    %fwrite(fileHandle,stf(i).numOfRays,'int');
    numOfBixels(i) = stf(i).numOfRays;

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
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

%% Call the OmpMC interface
%initialize waitbar
figureWait = waitbar(0,'OmpMC photon dose influence matrix calculation..');
fprintf('matRad: OmpMC photon dose calculation... ');

%run over all scenarios
for s = 1:dij.numOfScenarios
    fprintf('Scenario %d... ',s);
     
    ompMCgeo.isoCenter = [stf(:).isoCenter];
    
    
    
    dij.physicalDose{s} = matRad_ompInterface(cubeRho{s},cubeMatIx{s},ompMCgeo,ompMCsource,ompMCoptions);
end

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

end
