function [dij,ct,stf,pln,cst] = matRad_calcNeutronDoseMCNP(ct,stf,pln,cst,calcDoseDirect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad neutron dose calculation wrapper
%
% Neutron dose engine: Monte Carlo - MCNP6
% 
% call
%   dij = matRad_calcNeutronDoseMCNP(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: use predefined MLC field shape for calculation
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User’s Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 10/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prelude
matRad_cfg =  MatRad_Config.instance();

% initialize
matRad_calcDoseInit;

% Generate log-file
pathLog = strcat(matRad_cfg.matRadRoot,filesep, 'submodules', filesep, 'MCdoseEngineMCNP', filesep, 'logFile');

try diary(fullfile(pathLog, strcat(matRad_getTime4log, '_neutronDoseCalcuation')))
catch
    mkdir('logFile')
    diary(fullfile(pathLog, strcat(matRad_getTime4log, '_neutronDoseCalcuation')))
end

% Default: dose influence matrix computation
if ~exist('calcDoseDirect','var')
    calcDoseDirect = false;
    varHelper.calcDoseDirect = calcDoseDirect;
elseif exist('calcDoseDirect','var')
    if calcDoseDirect
        varHelper.calcDoseDirect = calcDoseDirect;
        stf(1).ray(1).energy = 'rssa';
    end
end


% Check if MCNP6 is installed
answer = questdlg('Is MCNP6 installed on your computer?', ...
    'MCNP6', ...
    'Yes', 'No', 'No');

switch answer
    case 'Yes'
        disp('*****')
        disp('Monte Carlo dose calculation enabled.')
        disp('*****')
        varHelper.runMCdoseCalc = 1;
    case 'No'
        disp('*****')
        disp(['MCNP runfiles will be created but no dose caculation will be performed.', newline, '(dij=zeros(size(ct.cubeHU)))'])
        disp('*****')
        varHelper.runMCdoseCalc = 0;
end
clear answer

% Load predefined conversion properties for tissue characterization
% according to CT values - elemental composition will be assigned according
% to these predefined tissue intervals later
disp('*****')
disp('Load pre-defined HU conversion properties and MCNP cross sections from conversionCT2tissue.mat.')
disp('Note: Modification of conversionCT2tissue.mat using generateVar_conversionCT2tissue.m.')
load('conversionCT2tissue.mat')
disp('*****')

%% Process CT data
% Check ct for MCNP Simulation
if ct.resolution.x ~= ct.resolution.y
    error('x- and y-resolution have to be equal for the simulation.')
end

% Set HU outside body to air, i.e. neglect everything outside body for the
% simulation
% Note: body structure is the only normal tissue structure that 
% has to be contured

pln.propMCNP.bodyStructureName = 'Body'; % Default name for body structure is 'Body'
try
    cstBodyIndex = matRad_findBodyStructureCST(cst, pln.propMCNP.bodyStructureName);
catch
    prompt = {'Please enter body structure name:'};
    dlgtitle = 'Find Body Structure';
    pln.propMCNP.bodyStructureName = inputdlg(prompt,dlgtitle);
    cstBodyIndex = matRad_findBodyStructureCST(cst, pln.propMCNP.bodyStructureName);
end

% Process HU values
disp('*****')
disp(['Properties from (scaled) HU loaded are: Minimum value: ',num2str(min(ct.cubeHU{1}, [], 'all')), ' and '])
disp(['Maximum value: ',num2str(max(ct.cubeHU{1}(cst{cstBodyIndex,4}{1}), [], 'all')), '.'])
if ~isfield(ct, 'dicomInfo') || ~isfield(ct.dicomInfo, 'RescaleSlope') || ~isfield(ct.dicomInfo, 'RescaleIntercept')
    matRad_cfg.dispWarning('No information on rescale slope and/or intercept provided in DICOM data. Calculation might crash...')
else
    disp(['Rescale slope from CT data is read to be: ',num2str(ct.dicomInfo.RescaleSlope), ' and rescale intercept is read to be ',num2str(ct.dicomInfo.RescaleIntercept),'.'])
end

disp('Please use question dialog to decide how to convert to scaled HU.')
disp('*****')

% Set values outside body to air
try
        cstTargetIndex = matRad_findTargetStructureCST(cst);
catch
    error('Target structure has to be set in matRad.')
end

maskNonBody = ones(ct.cubeDim);
bodyIdx = [cst{sort([cstTargetIndex, cstBodyIndex]),4}];
bodyIdx = unique(vertcat(bodyIdx{:}));

maskNonBody(bodyIdx) = 0;
ct.cube{1}(maskNonBody>0) = 0;
ct.cubeHU{1}(maskNonBody>0) = -1000;


% DICOM rescaling
answer = questdlg('Would you like to use DICOM rescale slope and intercept? If not, an offset of 1000 will be added to the HU values to get re-scale HUs.', ...
    'Use DICOM Info', ...
    'Yes', 'No', 'No');

switch answer
    case 'Yes'
        disp('*****')
        disp('You decided to use the following parameters to re-scale (scaled) HU data given in ct.cubeHU to HU in ct.cube.')
        disp(['Rescale HU: slope=', num2str(ct.dicomInfo.RescaleSlope), ' intercept=', num2str(ct.dicomInfo.RescaleIntercept)])
        disp('*****')
        ct.cubeHU{1} = ct.cubeHU{1}.*ct.dicomInfo.RescaleSlope - ct.dicomInfo.RescaleIntercept;
        %ct.cubeHU{1}(setxor(1:prod(ct.cubeDim), cst{cstBodyIndex,4}{1})) = 0;   % Set values outside body to air
        %ct.cubeHU{1}(ct.cubeHU{1}(:)<0) = 0;                                    % Clean up
    case 'No'
        disp('*****')
        disp('You decided not to use DICOM rescale slope and intercept.')
        disp('*****')
        ct.cubeHU{1} = ct.cubeHU{1} + 1000;
        %ct.cubeHU{1}(setxor(1:prod(ct.cubeDim), cst{cstBodyIndex,4}{1})) = 0;   % Set values outside body to air
        %ct.cubeHU{1}(ct.cubeHU{1}(:)<0) = 0;	                                % Clean up
end
clear answer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Information on Conversion                                                                        %
% Re-scale HU from intensity values given in DICOM by using 'rescale intercept' and 'rescale slope'           %
% Def. HU: HU = 1000*((mu - mu_water)/mu_water) -> HU_water = 0                                               %
% Note: For our purpose, HU values should start at zero.                                                      %         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Infere tissue characteristics from HU
disp('*****')
disp(['Properties from rescaled HU are: Minimum value: ',num2str(min(ct.cubeHU{1}, [], 'all')), ' and '])
disp(['Maximum value: ',num2str(max(ct.cubeHU{1}, [], 'all')), '.'])
disp('*****')
disp('Material types are assigned using the following HU intervals...')
maxHUbin_nonEmpty = 1;  % Find last non-empty entry
while ~isempty(binIntervals(maxHUbin_nonEmpty+1).HUbin) && (maxHUbin_nonEmpty <= size(binIntervals,2))
    maxHUbin_nonEmpty = maxHUbin_nonEmpty +1;
end
for i=1:maxHUbin_nonEmpty
    disp([binIntervals(i).name, ': ', num2str(binIntervals(i).HUbin(1)), ' to ', num2str(binIntervals(i).HUbin(2))]); 
end
disp('*****')

[cst, ct.tissueBin] = matRad_segmentationCTscan(ct.cubeHU{1}, ct.resolution, binIntervals, cst, cstBodyIndex, cstTargetIndex);

% Calculate density for CT voxels and resize afterwards -> caution: step-
% wise definition of conversion causes a difference concerning the ordering
% of the conversion and resize operation
disp('*****')
disp('Calculate density from CT data with density given in [g/cm^3]')
disp('*****')
ct.density{1} = hounsfield2density(ct.cubeHU{1})*1e-3; % rescale from kg/m^3 to g/cm^3
ct.density{1}(ct.density{1}<=(hounsfield2density(segVar.upperLimitAir)*1e-3)) = segVar.densityAir;

%% Rescale ct information to cm for MCNP runfile
varHelper.rescaleFactor = 1e-1; % conversion from mm to cm
ct.resolution.x_resized = ct.resolution.x*varHelper.rescaleFactor;
ct.resolution.y_resized = ct.resolution.y*varHelper.rescaleFactor;
ct.resolution.z_resized = ct.resolution.z*varHelper.rescaleFactor;


%% Create MCNP runfile blocks A and B
varHelper.simPropMCNP.loopCounter = false; % try to generate MCNP runfile with maximum 99999 elements for reasons of performance
varHelper.simPropMCNP.MCNP_limitNumberOfElements = 99999-1; % minus one since we need one integer for the source surface

[varHelper.simPropMCNP.control_makeTargetMCNP, varHelper.simPropMCNP.fileID_A, varHelper.simPropMCNP.fileID_B, varHelper.simPropMCNP.geometryOption] = ...
        matRad_makeTargetMCNP(ct, varHelper.simPropMCNP);

if ~varHelper.simPropMCNP.control_makeTargetMCNP
    varHelper.simPropMCNP.MCNP_limitNumberOfElements = 99999999-1;  % minus one since we need one integer for the source surface
    varHelper.simPropMCNP.loopCounter = true;
    [varHelper.simPropMCNP.control_makeTargetMCNP,~ ,~ ,~ ] = ...
        matRad_makeTargetMCNP(ct, varHelper.simPropMCNP);
end

if ~varHelper.simPropMCNP.control_makeTargetMCNP
    error('Number of defined elements for simulation is too high! MCNP6 only allows 99,999,999 cells and surfaces in total.')
end

%% Create MCNP runfile block C (source etc.)
% Set default number of particles 
defNPS = 1e6;
if isfield(pln, 'propMCNP') &&   isfield(pln.propMCNP, 'numberParticles')
    varHelper.simPropMCNP.numberParticles = pln.propMCNP.numberParticles;
else
    varHelper.simPropMCNP.numberParticles = defNPS;
    pln.propMCNP.numberParticles = defNPS;

end

% Get total number of bixels and write source card 
% Total number of bixel is counter for j in dij matrix

% C.1 Source
varHelper.totalNumberBixels = 0;
for counterField =1:size(stf,2)
    varHelper.simPropMCNP.counterField = counterField;
    for counterRay=1:stf(counterField).numOfRays
        varHelper.simPropMCNP.counterRay = counterRay;
        % Calculate position of MLC field in LPS system
        stf(counterField).ray(counterRay).rayPosMLC = stf(counterField).isoCenter + stf(counterField).ray(counterRay).rayPos + stf(counterField).sourcePoint;
        varHelper.totalNumberBixels = varHelper.totalNumberBixels + 1;
        
        % Generate source card for each bixel - see also description of
        % matRad^_makeSourceMCNP(...)
        [control_makeSourceMCNP, varHelper] = matRad_makeSourceMCNP(stf, pln, varHelper, counterField, counterRay);
    end
end

varHelper.simPropMCNP.sourceBlockNames = strings(1,varHelper.totalNumberBixels);
for cntList=1:varHelper.totalNumberBixels; varHelper.simPropMCNP.sourceBlockNames(cntList)=['blockC_source', int2str(cntList)]; end

% Open file to write rest of block C
% Source and the rest of block C input are seperated so the source
% positioning can be done easily w/o wasting time on redundant writing of
% the rest (like MODE and PHYS card) into a text file.

pathRunfiles = strcat(matRad_cfg.matRadRoot,filesep, 'submodules', filesep, 'MCdoseEngineMCNP', filesep, 'runfiles_tmp', filesep);
fileID_C_rest = fopen(strcat(pathRunfiles,'blockC_rest'), 'w');

% C.2 Physics and problem termination
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C.2: Physics\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');

matRad_definePhysicsMCNP(fileID_C_rest, pln, binIntervals, varHelper.simPropMCNP); % Define PHYS-card

% C.3 Materials
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C.3: Materials\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');

matIdentifierTxt = 'M%d $ %s \n';
matCompositionTxt = '           %d%s %.9f \n';

for counterMaterial = 1:size(binIntervals,2)
    fprintf(fileID_C_rest, matIdentifierTxt, counterMaterial, binIntervals(counterMaterial).name);
    fprintf(fileID_C_rest,'	   plib=14p\n');
    fprintf(fileID_C_rest,'	   hlib=70h\n');
    for counterComponent = 1:size(binIntervals(counterMaterial).ZAID,2)
            fprintf(fileID_C_rest,matCompositionTxt, binIntervals(counterMaterial).ZAID(counterComponent), crossSectionsLibrary(binIntervals(counterMaterial).crossSection(counterComponent),:), binIntervals(counterMaterial).percentageMass(counterComponent));
    end
end
clear counterMaterial; clear counterComponent;

% C.4 Tally
matRad_makeTallyMCNP(ct, pln, fileID_C_rest, binIntervals)

fclose(fileID_C_rest);

%% Concatenate all blocks to one runfile for each ray
matRad_concatenateRunfiles(varHelper, pathRunfiles);

%% Generate dij matrix
switch varHelper.runMCdoseCalc
    case 1
        dij = matRad_bixelDoseCalculatorMCNP(pathRunfiles, dij, stf, ct, pln, cst, binIntervals);
    case 0
        dij=zeros(size(ct.cubeHU));
end

%% Switch off diary
diary off
