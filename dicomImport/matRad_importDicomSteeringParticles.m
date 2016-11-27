function [stf, pln] = matRad_importDicomSteeringParticles(ct, pln, rtPlanFile)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to import a matRad stf struct from dicom RTPLAN data
% 
% call
%   [stf, pln] = matRad_importDicomSteeringParticles(ct, pln, rtPlanFile)
%
% input
%   ct:             ct imported by the matRad_importDicomCt function
%   cst:            cst created by the matRad_createCst function
%   rtPlanFile:   	name of RTPLAN DICOM file
%
% output
%   stf             matRad stf struct
%   pln:            matRad pln struct with meta information. Note that
%                   we also have here the pln output because pln.bixelWidth
%                   is determined here
%
% References
%   -
% Note
% not implemented - compensator. Fixed SAD.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load plan file
% load machine data

pln.machine = 'HIT';
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% RT Plan consists only on meta information
rtPlanInfo = dicominfo(rtPlanFile{1});
BeamSeq = rtPlanInfo.IonBeamSequence;
BeamSeqNames = fieldnames(BeamSeq);
% Number of Beams from plan
numOfBeamsPlan = length(pln.gantryAngles);

% use only the treatment beams
for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    try
        treatDelType = currBeamSeq.TreatmentDeliveryType;
        if ~strcmpi(treatDelType,'TREATMENT')
            BeamSeq = rmfield(BeamSeq,BeamSeqNames{i});
        end
    catch
        warning('Something went wrong while determining the type of the beam.');
    end
end

% reinitialize the BeamSeqNames and length, as the Seq itself is reduced.
BeamSeqNames = fieldnames(BeamSeq);

% remove empty ControlPointSequences
for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    ControlPointSeq      = currBeamSeq.IonControlPointSequence;
    ControlPointSeqNames = fieldnames(ControlPointSeq);
    numOfContrPointSeq = length(ControlPointSeqNames);
    for currContr = 1:numOfContrPointSeq
        currContrSeq = ControlPointSeq.(ControlPointSeqNames{currContr});
        if sum(currContrSeq.ScanSpotMetersetWeights) == 0
            ControlPointSeq = rmfield(ControlPointSeq,ControlPointSeqNames{currContr});
        end
    end
    BeamSeq.(BeamSeqNames{i}).IonControlPointSequence = ControlPointSeq;
end

% check if number of beams correspond
if ~isequal(length(BeamSeqNames),numOfBeamsPlan)
    warning('Number of beams from beamsequences do not correspond to number of Gantry Angles');
end

%% generate stf struct
% surfaceEntry = BeamSeq.Item_1.IonControlPointSequence.Item_1.SurfaceEntryPoint;

% Preallocate stf
stf(length(BeamSeqNames)).gantryAngle = [];
stf(length(BeamSeqNames)).couchAngle = [];
stf(length(BeamSeqNames)).bixelWidth = [];
stf(length(BeamSeqNames)).radiationMode = [];
stf(length(BeamSeqNames)).SAD = [];
stf(length(BeamSeqNames)).isoCenter = [];
stf(length(BeamSeqNames)).sourcePoint_bev = [];
stf(length(BeamSeqNames)).numOfRays = [];
stf(length(BeamSeqNames)).numOfBixelsPerRay = [];
stf(length(BeamSeqNames)).totalNumOfBixels = [];
stf(length(BeamSeqNames)).ray = [];

for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    ControlPointSeq      = currBeamSeq.IonControlPointSequence;
    stf(i).gantryAngle   = pln.gantryAngles(i);
    stf(i).couchAngle    = pln.couchAngles(i);
    stf(i).bixelWidth    = pln.bixelWidth;
    stf(i).radiationMode = pln.radiationMode;
    % there might be several SAD's, e.g. compensator?
    stf(i).SAD           = machine.meta.SAD;
    stf(i).isoCenter     = pln.isoCenter;
    stf(i).sourcePoint_bev = [0 -stf(i).SAD 0];
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient; use transpose matrices because we are
    % working with row vectors

    % Rotation around Z axis (gantry)
    rotMx_XY_T = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
                  -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
                                           0                         0 1];

    % Rotation around Y axis (couch)
    rotMx_XZ_T = [cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
                                         0 1                        0;
                  sind(pln.couchAngles(i)) 0  cosd(pln.couchAngles(i))];

    % Rotated Source point (1st gantry, 2nd couch)
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMx_XY_T*rotMx_XZ_T;

    % now loop over ControlPointSequences
    ControlPointSeqNames = fieldnames(ControlPointSeq);
    numOfContrPointSeq = length(ControlPointSeqNames);
    % create empty helper matrix
    StfTmp = zeros(0,5);
    for currContr = 1:numOfContrPointSeq
        currContrSeq = ControlPointSeq.(ControlPointSeqNames{currContr});
        % get energy, equal for all coming elements in the next loop
        currEnergy = currContrSeq.NominalBeamEnergy;
        % get focusValue
        currFocus = unique(currContrSeq.ScanningSpotSize);
        % get the Spotpositions
        numOfScanSpots = currContrSeq.NumberOfScanSpotPositions;
        % x is 1, 3, 5 ...; y 2, 4, 6,
        c1_help = currContrSeq.ScanSpotPositionMap(1:2:(2 * numOfScanSpots));
        c2_help = currContrSeq.ScanSpotPositionMap(2:2:(2 * numOfScanSpots));
        weight_help = currContrSeq.ScanSpotMetersetWeights;
        StfTmp = [StfTmp; c1_help c2_help (currEnergy * ones(numOfScanSpots,1)) weight_help (currFocus * ones(numOfScanSpots,1))];
    end
    
    % finds all unique rays and saves them in to the stf
    [RayPosTmp, ~, ic] = unique(StfTmp(:,1:2), 'rows');
    for j = 1:size(RayPosTmp,1)
        stf(i).ray(j).rayPos_bev = double([RayPosTmp(j,1) 0 RayPosTmp(j,2)]);
        stf(i).ray(j).energy = [];
        stf(i).ray(j).focusFWHM = [];
        stf(i).ray(j).focusIx = [];
        stf(i).ray(j).weight = [];
    end
    
    % saves all energies and weights to their corresponding ray
    for j = 1:size(StfTmp,1)
        k = ic(j);
        stf(i).ray(k).energy = [stf(i).ray(k).energy double(StfTmp(j,3))];
        stf(i).ray(k).focusFWHM = [stf(i).ray(k).focusFWHM double(StfTmp(j,5))];
        stf(i).ray(k).weight = [stf(i).ray(k).weight double(StfTmp(j,4)) / 1e6];
    end
    
    
    % getting some information of the rays
    % clean up energies, so they appear only one time per energy
    numOfRays = size(stf(i).ray,2);
    for l = 1:numOfRays
        stf(i).ray(l).energy = unique(stf(i).ray(l).energy);
        stf(i).ray(l).targetPoint_bev = [2*stf(i).ray(l).rayPos_bev(1) ...
                                         machine.meta.SAD ...
                                         2*stf(i).ray(l).rayPos_bev(3)];
    end
    stf(i).numOfRays = numel(stf(i).ray);  
    
    % save total number of bixels
    numOfBixels = 0;
    for j = 1:numel(stf(i).ray)
        numOfBixels = numOfBixels + numel(stf(i).ray(j).energy);
        stf(i).numOfBixelsPerRay(j) = numel(stf(i).ray(j).energy);
%         w = [w stf(currBeam).ray(j).weight];
    end
    
    stf(i).totalNumOfBixels = numOfBixels;
    
    % get bixelwidth
    bixelWidth_help = zeros(size(stf(i).ray,2),2);
    for j = 1:stf(i).numOfRays
        bixelWidth_help(j,1) = stf(i).ray(j).rayPos_bev(1);
        bixelWidth_help(j,2) = stf(i).ray(j).rayPos_bev(3);        
    end
    bixelWidth_help1 = unique(round(1e3*bixelWidth_help(:,1))/1e3,'sorted');
    bixelWidth_help2 = unique(round(1e3*bixelWidth_help(:,2))/1e3,'sorted');
    
    bixelWidth = unique([unique(diff(bixelWidth_help1)) unique(diff(bixelWidth_help2))]);
    
    if numel(bixelWidth) == 1
        stf(i).bixelWidth = bixelWidth;
    else
        stf(i).bixelWidth = NaN;
    end
    
    % compute coordinates in lps coordinate system, i.e. rotate beam
    % geometry around fixed patient
    
    % Rotation around Z axis (gantry)
    rotMx_XY_rotated = [ cosd(pln.gantryAngles(i)) sind(pln.gantryAngles(i)) 0;
        -sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
        0                         0 1];
    
    % Rotation around Y axis (couch)
    rotMx_XZ_rotated = [ cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
        0 1                        0;
        sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % Rotated Source point, first needs to be rotated around gantry, and then
    % couch.
    stf(i).sourcePoint =  stf(i).sourcePoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMx_XY_rotated*rotMx_XZ_rotated;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMx_XY_rotated*rotMx_XZ_rotated;   
        stf(i).ray(j).SSD         = NaN;
    end
    
    % SSD
    DensityThresholdSSD = 0.05;
    for j = 1:stf(i).numOfRays
        [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                             ct.resolution, ...
                             stf(i).sourcePoint, ...
                             stf(i).ray(j).targetPoint, ...
                             ct.cube);

            ixSSD = find(rho{1} > DensityThresholdSSD,1,'first');

            if isempty(ixSSD)== 1
                warning('Surface for SSD calculation starts directly in first voxel of CT\n');
            end
            
            % calculate SSD
            stf(i).ray(j).SSD = double(2 * stf(i).SAD * alpha(ixSSD));
            % book keeping & calculate focus index
            stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
    end
    
    % use the original machine energies
    for j = 1:stf(i).numOfRays
        % loop over all energies
        numOfEnergy = length(stf(i).ray(j).energy);
        for k = 1:numOfEnergy
            energyIndex = find(abs([machine.data(:).energy]-stf(i).ray(j).energy(k))<10^-3);
            if ~isempty(energyIndex)
                stf(i).ray(j).energy(k) = machine.data(energyIndex).energy;
            else
                error('No match between imported and machine data. Maybe wrong machine loaded.');
            end
        end
    end
    
    % get focusIx instead of focusFWHM
    for j = 1:stf(i).numOfRays
        % loop over all energies
        numOfEnergy = length(stf(i).ray(j).energy);
        for k = 1:numOfEnergy
            energyTemp = stf(i).ray(j).energy(k);
            focusFWHM = stf(i).ray(j).focusFWHM(k);
            energyIxTemp = find([machine.data.energy] == energyTemp);
            focusIxTemp = find(abs([machine.data(energyIxTemp).initFocus.SisFWHMAtIso] - focusFWHM )< 10^-3);
            stf(i).ray(j).focusIx(k) = focusIxTemp;
            stf(i).ray(j).focusFWHM(k) = machine.data(energyIxTemp).initFocus.SisFWHMAtIso(stf(i).ray(j).focusIx(k));
        end
    end
    
    stf(i).timeStamp = datestr(clock);
    
end

if any(isnan([stf(:).bixelWidth])) || numel(unique([stf(:).bixelWidth])) > 1
    pln.bixelWidth = NaN;
else
    pln.bixelWidth = stf(1).bixelWidth;
end

end
