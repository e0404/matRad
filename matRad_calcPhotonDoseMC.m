function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,nCasePerBixel,visBool)
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
%   -
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


matRad_cfg =  MatRad_Config.instance();

tic

% disable visualiazation by default
if nargin < 6
    visBool = false;
end

if nargin < 5
    nCasePerBixel = matRad_cfg.propMC.ompMC_defaultHistories;
    matRad_cfg.dispInfo('Using default number of Histories per Bixel: %d\n',nCasePerBixel);
end

fileFolder = fileparts(mfilename('fullpath'));

if ~matRad_checkMexFileExists('omc_matrad') %exist('matRad_ompInterface','file') ~= 3    
    matRad_cfg.dispWarning('Compiled mex interface not found. Trying to compile the ompMC interface on the fly!');    
    try        
        matRad_compileOmpMCInterface();
    catch MException       
        matRad_cfg.dispError('Could not find/generate mex interface for MC dose calculation.\nCause of error:\n%s\n Please compile it yourself (preferably with OpenMP support).',MException.message);
    end
end

matRad_calcDoseInit;

% gaussian filter to model penumbra from (measured) machine output / see diploma thesis siggel 4.1.2
if isfield(machine.data,'penumbraFWHMatIso')
    penumbraFWHM = machine.data.penumbraFWHMatIso;
else
    penumbraFWHM = 5;
    matRad_cfg.dispWarning('photon machine file does not contain measured penumbra width in machine.data.penumbraFWHMatIso. Assuming 5 mm.');
end

sourceFWHM = penumbraFWHM * machine.meta.SCD/(machine.meta.SAD - machine.meta.SCD);
sigmaGauss = sourceFWHM / sqrt(8*log(2)); % [mm] 

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfBixels,1);
dij.rayNum   = NaN*ones(dij.totalNumOfBixels,1);
dij.beamNum  = NaN*ones(dij.totalNumOfBixels,1);

dij.numHistoriesPerBeamlet = nCasePerBixel;

omcFolder = [matRad_cfg.matRadRoot filesep 'ompMC'];
%omcFolder = [matRad_cfg.matRadRoot filesep 'submodules' filesep 'ompMC'];

%% Setup OmpMC options / parameters

%display options
ompMCoptions.verbose = matRad_cfg.logLevel - 1;

% start MC control          
ompMCoptions.nHistories = nCasePerBixel;
ompMCoptions.nSplit = 20;
ompMCoptions.nBatches = 10;
ompMCoptions.randomSeeds = [97 33];

%start source definition      
ompMCoptions.spectrumFile       = [omcFolder filesep 'spectra' filesep 'mohan6.spectrum'];
ompMCoptions.monoEnergy         = 0.1; 
ompMCoptions.charge             = 0;
ompMCoptions.sourceGeometry     = 'gaussian';
ompMCoptions.sourceGaussianWidth = 0.1*sigmaGauss;
                                                                    
% start MC transport
ompMCoptions.dataFolder   = [omcFolder filesep 'data' filesep];
ompMCoptions.pegsFile     = [omcFolder filesep 'pegs4' filesep '700icru.pegs4dat'];
ompMCoptions.pgs4formFile = [omcFolder filesep 'pegs4' filesep 'pgs4form.dat'];

ompMCoptions.global_ecut = 0.7;
ompMCoptions.global_pcut = 0.010; 

% Relative Threshold for dose
ompMCoptions.relDoseThreshold = 1 - matRad_cfg.propDoseCalc.defaultLateralCutOff;

% Output folders
ompMCoptions.outputFolder = [omcFolder filesep 'output' filesep];

% Create Material Density Cube
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

    HUcube{s} =  matRad_interp3(dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,ct.cubeHU{s}, ...
                            dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
    
    % projecting out of bounds HU values where necessary
    if max(HUcube{s}(:)) > material{end,3}
        matRad_cfg.dispWarning('Projecting out of range HU values');
        HUcube{s}(HUcube{s}(:) > material{end,3}) = material{end,3};
    end
    if min(HUcube{s}(:)) < material{1,2}
        matRad_cfg.dispWarning('Projecting out of range HU values');
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

scale = 10; % to convert to cm

ompMCgeo.xBounds = (dij.doseGrid.resolution.y * (0.5 + [0:dij.doseGrid.dimensions(1)])) ./ scale;
ompMCgeo.yBounds = (dij.doseGrid.resolution.x * (0.5 + [0:dij.doseGrid.dimensions(2)])) ./ scale;
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
useCornersSCD = true; %false -> use ISO corners

numOfBixels = [stf(:).numOfRays];
beamSource = zeros(dij.numOfBeams, 3);

bixelCorner = zeros(dij.totalNumOfBixels,3);
bixelSide1 = zeros(dij.totalNumOfBixels,3);
bixelSide2 = zeros(dij.totalNumOfBixels,3);

counter = 0;

for i = 1:dij.numOfBeams % loop over all beams
   
    % define beam source in physical coordinate system in cm
    beamSource(i,:) = (stf(i).sourcePoint + stf(i).isoCenter)/10;

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        counter = counter + 1;
        
        dij.beamNum(counter)  = i;
        dij.rayNum(counter)   = j;
        dij.bixelNum(counter) = j;
        
        if useCornersSCD
            beamletCorners = stf(i).ray(j).rayCorners_SCD;
        else    
            beamletCorners = stf(i).ray(j).beamletCornersAtIso;
        end
        
        % get bixel corner and delimiting vectors.
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        currCorner = (beamletCorners(1,:) + stf(i).isoCenter) ./ scale;
        bixelCorner(counter,:) = currCorner;
        bixelSide1(counter,:) = (beamletCorners(2,:) + stf(i).isoCenter) ./ scale - currCorner;
        bixelSide2(counter,:) = (beamletCorners(4,:) + stf(i).isoCenter) ./ scale - currCorner;
        
        if visBool
            for k = 1:4
                currCornerVis = (beamletCorners(k,:) + stf(i).isoCenter)/10;
                % rays connecting source and ray corner
                plot3([beamSource(i,1) currCornerVis(1)],[beamSource(i,2) currCornerVis(2)],[beamSource(i,3) currCornerVis(3)],'b')
                % connection between corners
                lRayCorner = (beamletCorners(mod(k,4) + 1,:) + stf(i).isoCenter)/10;
                plot3([lRayCorner(1) currCornerVis(1)],[lRayCorner(2) currCornerVis(2)],[lRayCorner(3) currCornerVis(3)],'r')
            end
        end
        
    end
        
end

ompMCsource.nBeams = dij.numOfBeams;
ompMCsource.iBeam = dij.beamNum(:);

% Switch x and y directions to match ompMC cs.
ompMCsource.xSource = beamSource(:,2);
ompMCsource.ySource = beamSource(:,1);
ompMCsource.zSource = beamSource(:,3);

ompMCsource.nBixels = sum(numOfBixels(:));
ompMCsource.xCorner = bixelCorner(:,2);
ompMCsource.yCorner = bixelCorner(:,1);
ompMCsource.zCorner = bixelCorner(:,3);

ompMCsource.xSide1 = bixelSide1(:,2);
ompMCsource.ySide1 = bixelSide1(:,1);
ompMCsource.zSide1 = bixelSide1(:,3);

ompMCsource.xSide2 = bixelSide2(:,2);
ompMCsource.ySide2 = bixelSide2(:,1);
ompMCsource.zSide2 = bixelSide2(:,3);

if visBool
    plot3(ompMCsource.ySource,ompMCsource.xSource,ompMCsource.zSource,'rx')
end

%% Call the OmpMC interface

%ompMC for matRad returns dose/history * nHistories. 
% This factor calibrates to 1 Gy in a %(5x5)cm^2 open field (1 bixel) at 
% 5cm depth for SSD = 900 which corresponds to the calibration for the 
% analytical base data.
absCalibrationFactor = 3.49056 * 1e12; %Approximate!

%Now we have to calibrate to the the beamlet width.
absCalibrationFactor = absCalibrationFactor * (pln.propStf.bixelWidth/50)^2;

matRad_cfg.dispInfo('matRad: OmpMC photon dose calculation... \n');

outputVariance = matRad_cfg.propMC.ompMC_defaultOutputVariance;

if isfield(pln,'propMC') && isfield(pln.propMC,'outputVariance')
    outputVariance = pln.propMC.outputVariance;
end


%run over all scenarios
for s = 1:dij.numOfScenarios
    ompMCgeo.isoCenter = [stf(:).isoCenter];
    
    %Run the Monte Carlo simulation and catch  possible mex-interface
    %issues
    try
        %If we ask for variance, a field in the dij will be filled
        if outputVariance
            [dij.physicalDose{s},dij.physicalDose_MCvar{s}] = omc_matrad(cubeRho{s},cubeMatIx{s},ompMCgeo,ompMCsource,ompMCoptions);
        else
            [dij.physicalDose{s}] = omc_matrad(cubeRho{s},cubeMatIx{s},ompMCgeo,ompMCsource,ompMCoptions);
        end
    catch ME
        errorString = [ME.message '\nThis error was thrown by the MEX-interface of ompMC.\nMex interfaces can raise compatability issues which may be resolved by compiling them by hand directly on your particular system.'];
        matRad_cfg.dispError(errorString);
    end
    
    %Calibrate the dose with above factor
    dij.physicalDose{s} = dij.physicalDose{s} * absCalibrationFactor;
    if isfield(dij,'physicalDose_MCvar')
        dij.physicalDose_MCvar{s} = dij.physicalDose_MCvar{s} * absCalibrationFactor^2;
    end
end

matRad_cfg.dispInfo('matRad: MC photon dose calculation done!\n');
matRad_cfg.dispInfo(evalc('toc'));

try
    % wait 0.1s for closing all waitbars
    allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar');
    delete(allWaitBarFigures);
    pause(0.1);
catch
end

end
