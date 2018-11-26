function dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,nCasePerBixel,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad vmc++ photon dose calculation wrapper
% 
% call
%   dij = matRad_calcPhotonDoseMC(ct,stf,pln,cst,nCasePerBixel)
%
% input
%   ct:                         matRad ct struct
%   stf:                        matRad steering information struct
%   pln:                        matRad plan meta information struct
%   cst:                        matRad cst struct
%   nCasePerBixel:              number of photons simulated per bixel
%   calcDoseDirect:             boolian switch to bypass dose influence matrix
%                               computation and directly calculate dose; only makes
%                               sense in combination with matRad_calcDoseDirect.m%
% output
%   dij:                        matRad dij struct
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
end

if ~isdeployed % only if _not_ running as standalone    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'vmc++'))
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
if calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

bixelNum = NaN*ones(dij.totalNumOfBixels,1);
rayNum   = NaN*ones(dij.totalNumOfBixels,1);
beamNum  = NaN*ones(dij.totalNumOfBixels,1);

% Allocate space for dij.physicalDose sparse matrix
for i = 1:dij.numOfScenarios
    dij.physicalDose{i} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
end

% set consistent random seed (enables reproducibility)
rng(0);

% set number of photons simulated per bixel and number of parallel MC simulations if not specified by user
if nargin < 5
    nCasePerBixel              = 5000;
end
    
% set relative dose cutoff for storage in dose influence matrix
relDoseCutoff = 10^(-3);

% set absolute calibration factor
% CALCULATION
% absolute_calibration_factor = 1/D(depth = 100,5mm) -> D(depth = 100,5mm) = 1Gy
% SETUP
% SAD = 1000mm, SCD = 500mm, bixelWidth = 5mm, IC = [240mm,240mm,240mm]
% fieldsize@IC = 105mm x 105mm, phantomsize = 81 x 81 x 81 = 243mm x 243mm x 243mm
% rel_Dose_cutoff = 10^(-3), ncase = 500000/bixel
absCalibrationFactorVmc = 99.818252282632300;

% set general vmc++ parameters
% 1 source
mcOptions.beamletSource.myName       = 'source 1';                        % name of source
mcOptions.beamletSource.monitorUnits = 1;                                 
mcOptions.beamletSource.spectrum     = ['./spectra/var_6MV.spectrum'];    % energy spectrum source (only used if no mono-Energy given)
mcOptions.beamletSource.charge       = 0;                                 % charge (-1,0,1)
% 2 transport parameter
mcOptions.McParameter.automatic_parameter = 'yes';                        % if yes, automatic transport parameters are used
% 3 MC control
mcOptions.McControl.ncase  = nCasePerBixel;                               % number of histories
mcOptions.McControl.nbatch = 10;                                          % number of batches
% 4 variance reduction
mcOptions.varianceReduction.repeatHistory      = 0.251;
mcOptions.varianceReduction.splitPhotons       = 'yes';   
mcOptions.varianceReduction.photonSplitFactor = -40;  
% 5 quasi random numbers
mcOptions.quasi.base      = 2;                                                 
mcOptions.quasi.dimension = 60;                                             
mcOptions.quasi.skip      = 1;                                              
% 6 geometry
mcOptions.geometry.XyzGeometry.methodOfInput = 'CT-PHANTOM';              % input method ('CT-PHANTOM', 'individual', 'groups') 
mcOptions.geometry.XyzGeometry.Ct            = 'CT';                      % name of geometry
mcOptions.geometry.XyzGeometry.CtFile        = ['./phantoms/matRad_CT.ct']; % path of density matrix (only needed if input method is 'CT-PHANTOM')
% 7 scoring manager
mcOptions.scoringOptions.startInGeometry               = 'CT';            % geometry in which partciles start their transport
mcOptions.scoringOptions.doseOptions.scoreInGeometries = 'CT';            % geometry in which dose is recorded
mcOptions.scoringOptions.doseOptions.scoreDoseToWater  = 'yes';           % if yes output is dose to water
mcOptions.scoringOptions.outputOptions.name            = 'CT';            % geometry for which dose output is created (geometry has to be scored)
mcOptions.scoringOptions.outputOptions.dumpDose        = 2;               % output format (1: format=float, Dose + deltaDose; 2: format=short int, Dose)

% export CT cube as binary file for vmc++
%matRad_exportCtVmc(ct, fullfile(phantomPath, 'matRad_CT.ct'));

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

writeCounter                  = 0;

% open beamlet file for reading
fileHandle = fopen('beamlets.dat','w');
fwrite(fileHandle,dij.numOfBeams,'int');

% visualization
clf
hold on
axis equal
% ct box
ctCorner1 = [ct.resolution.x ct.resolution.y ct.resolution.z]/2/10;
ctCorner2 = ([ct.cubeDim + .5] .* [ct.resolution.x ct.resolution.y ct.resolution.z])/10;
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

plot3(pln.propStf.isoCenter(1)/10,pln.propStf.isoCenter(2)/10,pln.propStf.isoCenter(3)/10,'gx')

xlabel('x [cm]')
ylabel('y [cm]')
zlabel('z [cm]')

rotate3d on

fprintf('matRad: MC photon dose calculation... ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams % loop over all beams
       
   % remember beam and bixel number
    if calcDoseDirect
        dij.beamNum(i)    = i;
        dij.rayNum(i)     = i;
        dij.bixelNum(i)   = i;
    end

    % define beam source in physical coordinate system in cm
    beamSource = (stf(i).sourcePoint + stf(i).isoCenter)/10;
    fwrite(fileHandle,beamSource,'double');   
    
    % write number of beamlets into file
    fwrite(fileHandle,stf(i).numOfRays,'int');

    for j = 1:stf(i).numOfRays % loop over all rays / for photons we only have one bixel per ray!
        
        writeCounter = writeCounter + 1;

        % create different seeds for every bixel
        mcOptions.McControl.rngSeeds = [randi(30000),randi(30000)];

        % remember beam and bixel number
        if ~calcDoseDirect
           dij.beamNum(writeCounter)  = i;
           dij.rayNum(writeCounter)   = j;
           dij.bixelNum(writeCounter) = j;
        end
        beamNum(writeCounter)  = i;
        rayNum(writeCounter)   = j;
        bixelNum(writeCounter) = j;
        
        % set ray specific vmc++ parameters
        % a) change coordinate system (Isocenter cs-> physical cs) and units mm -> cm
        for k = 1:4
            currRayCorner = (stf(i).ray(j).beamletCornersAtIso(k,:) + stf(i).isoCenter)/10;
            fwrite(fileHandle,currRayCorner,'double');   
            
            % rays connecting source and ray corner
            plot3([beamSource(1) currRayCorner(1)],[beamSource(2) currRayCorner(2)],[beamSource(3) currRayCorner(3)],'y')
            % connection between corners
            l = mod(k,4) + 1;
            lRayCorner = (stf(i).ray(j).beamletCornersAtIso(l,:) + stf(i).isoCenter)/10;
            plot3([lRayCorner(1) currRayCorner(1)],[lRayCorner(2) currRayCorner(2)],[lRayCorner(3) currRayCorner(3)],'r')
            
        end
    end
end

fclose(fileHandle);

% delete temporary files

try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
