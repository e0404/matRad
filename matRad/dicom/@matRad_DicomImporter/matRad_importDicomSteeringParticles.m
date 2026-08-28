function obj = matRad_importDicomSteeringParticles(obj)
% matRad function to import a matRad stf struct from dicom RTPLAN data
%
% In your object, there must be properties that contain:
%   - ct imported by the matRad_importDicomCt function;
%   - matRad pln structure with meta information;
%   - name of RTPLAN DICOM file.
%
% Output - matRad stf and pln structures.
% Note: pln is input and output since pln.bixelWidth is determined here.
%
% call:
%   obj = matRad_importDicomSteeringParticles(obj)
%
%
% References
%   -
%
% Note
% not implemented - compensator. Fixed SAD.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

%% load plan file
% load machine data

matRad_cfg = MatRad_Config.instance();
matRad_checkEnvDicomRequirements(matRad_cfg.env);

dlgBaseDataText = ['Import steering information from DICOM Plan.', 'Choose corresponding matRad base data for ', ...
                   obj.pln.radiationMode, '.'];
% messagebox only necessary for non windows users
if ~ispc
    uiwait(helpdlg(dlgBaseDataText, ['DICOM import - ', obj.pln.radiationMode, ' base data']));
end
[fileName, pathName] = uigetfile([matRad_cfg.matRadSrcRoot filesep 'basedata' filesep '*.mat'], dlgBaseDataText);
load([pathName filesep fileName]);

ix = find(fileName == '_');
obj.pln.machine = fileName(ix(1) + 1:end - 4);

% RT Plan consists only on meta information
if matRad_cfg.isOctave || verLessThan('matlab', '9')
    rtPlanInfo = dicominfo(obj.importFiles.rtplan{1});
else
    rtPlanInfo = dicominfo(obj.importFiles.rtplan{1}, 'UseDictionaryVR', true);
end
BeamSeq = rtPlanInfo.IonBeamSequence;
BeamSeqNames = fieldnames(BeamSeq);
% Number of Beams from plan
numOfBeamsPlan = length(obj.pln.propStf.gantryAngles);

% use only the treatment beams
for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    try
        treatDelType = currBeamSeq.TreatmentDeliveryType;
        if ~strcmpi(treatDelType, 'TREATMENT')
            BeamSeq = rmfield(BeamSeq, BeamSeqNames{i});
        end
    catch
        warning('Something went wrong while determining the type of the beam.');
    end
end

% reinitialize the BeamSeqNames and length, as the Seq itself is reduced.
BeamSeqNames = fieldnames(BeamSeq);

% remove empty ControlPointSequences (PBS beams only)
for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    if strcmp(currBeamSeq.ScanMode, 'NONE')
        % Passive scattering beam - no scan spot pruning applicable
        continue
    end
    ControlPointSeq      = currBeamSeq.IonControlPointSequence;
    ControlPointSeqNames = fieldnames(ControlPointSeq);
    numOfContrPointSeq = length(ControlPointSeqNames);
    for currContr = 1:numOfContrPointSeq
        currContrSeq = ControlPointSeq.(ControlPointSeqNames{currContr});
        if sum(currContrSeq.ScanSpotMetersetWeights) == 0
            ControlPointSeq = rmfield(ControlPointSeq, ControlPointSeqNames{currContr});
        end
    end
    BeamSeq.(BeamSeqNames{i}).IonControlPointSequence = ControlPointSeq;
end

% check if number of beams correspond
if ~isequal(length(BeamSeqNames), numOfBeamsPlan)
    warning('Number of beams from beamsequences do not correspond to number of Gantry Angles');
end

%% generate stf struct
% surfaceEntry = BeamSeq.Item_1.IonControlPointSequence.Item_1.SurfaceEntryPoint;

% Preallocate stf
obj.stf(length(BeamSeqNames)).gantryAngle = [];
obj.stf(length(BeamSeqNames)).couchAngle = [];
obj.stf(length(BeamSeqNames)).bixelWidth = [];
obj.stf(length(BeamSeqNames)).radiationMode = [];
obj.stf(length(BeamSeqNames)).SAD = [];
obj.stf(length(BeamSeqNames)).isoCenter = [];
obj.stf(length(BeamSeqNames)).sourcePoint_bev = [];
obj.stf(length(BeamSeqNames)).numOfRays = [];
obj.stf(length(BeamSeqNames)).numOfBixelsPerRay = [];
obj.stf(length(BeamSeqNames)).totalNumOfBixels = [];
obj.stf(length(BeamSeqNames)).ray = [];
obj.stf(length(BeamSeqNames)).deliveryType = [];

for i = 1:length(BeamSeqNames)
    currBeamSeq = BeamSeq.(BeamSeqNames{i});
    ControlPointSeq      = currBeamSeq.IonControlPointSequence;
    obj.stf(i).gantryAngle   = obj.pln.propStf.gantryAngles(i);
    obj.stf(i).couchAngle    = obj.pln.propStf.couchAngles(i);
    obj.stf(i).bixelWidth    = obj.pln.propStf.bixelWidth;
    obj.stf(i).radiationMode = obj.pln.radiationMode;
    % SAD - PBS beams provide two virtual source distances, passive beams a single scalar
    if isfield(currBeamSeq, 'VirtualSourceAxisDistances')
        obj.stf(i).SAD_x = currBeamSeq.VirtualSourceAxisDistances(1);
        obj.stf(i).SAD_y = currBeamSeq.VirtualSourceAxisDistances(2);
    elseif isfield(currBeamSeq, 'SourceAxisDistance')
        obj.stf(i).SAD_x = currBeamSeq.SourceAxisDistance;
        obj.stf(i).SAD_y = currBeamSeq.SourceAxisDistance;
    else
        matRad_cfg.dispWarning('Could not determine SAD for beam %d. Setting to NaN.', i);
        obj.stf(i).SAD_x = NaN;
        obj.stf(i).SAD_y = NaN;
    end
    % stf(i).SAD           = machine.meta.SAD; %we write the SAD later when we check machine match
    % stf(i).sourcePoint_bev = [0 -stf(i).SAD 0];
    obj.stf(i).isoCenter     = obj.pln.propStf.isoCenter(i, :);

    if strcmp(currBeamSeq.ScanMode, 'NONE')
        % Passive scattering beam: import geometry and represent as a single central ray/bixel.
        % The beam is not yet calculable in matRad, but the stf structure is kept consistent
        % so it can be integrated with future compensator-based pencil-beam adaptations.
        matRad_cfg.dispWarning('Beam %d uses passive scattering (ScanMode NONE). Imported as single central ray, but not yet calculable in matRad!', i);
        obj.stf(i).deliveryType = 'passiveScattering';

        % Extract energy, total weight and range shifter from control points
        ControlPointSeqNames = fieldnames(ControlPointSeq);
        beamEnergy   = NaN;
        beamWeight   = 1;
        rashiID              = 0;
        rashiWeThickness     = 0;
        rashiIsoRangeDist    = 0;
        for currContr = 1:length(ControlPointSeqNames)
            currContrSeq = ControlPointSeq.(ControlPointSeqNames{currContr});
            if isnan(beamEnergy) && isfield(currContrSeq, 'NominalBeamEnergy')
                beamEnergy = double(currContrSeq.NominalBeamEnergy);
            end
            if isfield(currContrSeq, 'FinalCumulativeMetersetWeight')
                beamWeight = double(currContrSeq.FinalCumulativeMetersetWeight);
            end
            if isfield(currContrSeq, 'RangeShifterSettingsSequence')
                rashiID           = currContrSeq.RangeShifterSettingsSequence.Item_1.ReferencedRangeShifterNumber;
                rashiWeThickness  = currContrSeq.RangeShifterSettingsSequence.Item_1.RangeShifterWaterEquivalentThickness;
                rashiIsoRangeDist = currContrSeq.RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance;
            end
        end

        % Single central ray at isocenter (BEV coordinates)
        obj.stf(i).ray(1).rayPos_bev  = [0 0 0];
        obj.stf(i).ray(1).energy      = beamEnergy;
        obj.stf(i).ray(1).focusFWHM   = NaN;   % not applicable for passive delivery
        obj.stf(i).ray(1).focusIx     = 1;
        obj.stf(i).ray(1).weight      = beamWeight / 1e6;
        obj.stf(i).ray(1).rangeShifter(1).ID                  = rashiID;
        obj.stf(i).ray(1).rangeShifter(1).eqThickness         = rashiWeThickness;
        obj.stf(i).ray(1).rangeShifter(1).sourceRashiDistance = obj.stf(i).SAD - rashiIsoRangeDist;

        obj.stf(i).numOfRays         = 1;
        obj.stf(i).numOfBixelsPerRay = 1;
        obj.stf(i).totalNumOfBixels  = 1;
        obj.stf(i).bixelWidth        = NaN;
    else
        obj.stf(i).deliveryType = 'scanned';

        % now loop over ControlPointSequences
        ControlPointSeqNames = fieldnames(ControlPointSeq);
        numOfContrPointSeq = length(ControlPointSeqNames);
        % create empty helper matrix
        temporarySteering = zeros(0, 8);
        for currContr = 1:numOfContrPointSeq
            currContrSeq = ControlPointSeq.(ControlPointSeqNames{currContr});
            % get energy, equal for all coming elements in the next loop
            if ~isfield(currContrSeq, 'NominalBeamEnergy')
                continue
            end
            currEnergy = currContrSeq.NominalBeamEnergy;
            % get focusValue and spot positions
            currFocus = unique(currContrSeq.ScanningSpotSize);
            numOfScanSpots = currContrSeq.NumberOfScanSpotPositions;
            % x is 1, 3, 5 ...; y 2, 4, 6,
            c1_help = currContrSeq.ScanSpotPositionMap(1:2:(2 * numOfScanSpots));
            c2_help = currContrSeq.ScanSpotPositionMap(2:2:(2 * numOfScanSpots));
            weight_help = currContrSeq.ScanSpotMetersetWeights;

            if isfield(currContrSeq, 'RangeShifterSettingsSequence')
                % rangeshifter identification
                rashiID = currContrSeq.RangeShifterSettingsSequence.Item_1.ReferencedRangeShifterNumber;
                % rangeshifter waterequivalent thickness
                rashiWeThickness = currContrSeq.RangeShifterSettingsSequence.Item_1.RangeShifterWaterEquivalentThickness;
                % rangeshifter isocenter to range shifter distance
                rashiIsoRangeDist = currContrSeq.RangeShifterSettingsSequence.Item_1.IsocenterToRangeShifterDistance;
            elseif currContr == 1
                rashiID = 0;
                rashiWeThickness = 0;
                rashiIsoRangeDist = 0;
            else
                % in this case range shifter settings has not changed between this
                % and previous control sequence, so reuse values.
            end
            temporarySteering = [temporarySteering; c1_help c2_help ...
                                 (currEnergy * ones(numOfScanSpots, 1)) weight_help (currFocus * ones(numOfScanSpots, 1)) ...
                                 (rashiID * ones(numOfScanSpots, 1)) (rashiWeThickness * ones(numOfScanSpots, 1)) (rashiIsoRangeDist * ones(numOfScanSpots, 1))];
        end

        % finds all unique rays and saves them in to the stf
        [RayPosTmp, ~, ic] = unique(temporarySteering(:, 1:2), 'rows');
        clear ray;
        for j = 1:size(RayPosTmp, 1)
            obj.stf(i).ray(j).rayPos_bev = double([RayPosTmp(j, 1) 0 RayPosTmp(j, 2)]);
            obj.stf(i).ray(j).energy = [];
            obj.stf(i).ray(j).focusFWHM = [];
            obj.stf(i).ray(j).focusIx = [];
            obj.stf(i).ray(j).weight = [];
            obj.stf(i).ray(j).rangeShifter = struct();
            ray(j).ID = [];
            ray(j).eqThickness = [];
            ray(j).sourceRashiDistance = [];
        end

        % saves all energies and weights to their corresponding ray
        for j = 1:size(temporarySteering, 1)
            k = ic(j);
            obj.stf(i).ray(k).energy = [obj.stf(i).ray(k).energy double(temporarySteering(j, 3))];
            obj.stf(i).ray(k).focusFWHM = [obj.stf(i).ray(k).focusFWHM double(temporarySteering(j, 5))];
            obj.stf(i).ray(k).weight = [obj.stf(i).ray(k).weight double(temporarySteering(j, 4)) / 1e6];
            % helpers to construct something like a(:).b = c.b(:) after this loop
            ray(k).ID = [ray(k).ID double(temporarySteering(j, 6))];
            ray(k).eqThickness = [ray(k).eqThickness double(temporarySteering(j, 7))];
            ray(k).sourceRashiDistance = [ray(k).sourceRashiDistance double(temporarySteering(j, 8))];
        end

        % reassign to preserve data structure
        for j = 1:numel(ray)
            for k = 1:numel(ray(j).ID)
                obj.stf(i).ray(j).rangeShifter(k).ID = ray(j).ID(k);
                obj.stf(i).ray(j).rangeShifter(k).eqThickness = ray(j).eqThickness(k);
                obj.stf(i).ray(j).rangeShifter(k).sourceRashiDistance = obj.stf(i).SAD - ray(j).sourceRashiDistance(k);
            end
        end

        % getting some information of the rays
        % clean up energies, so they appear only one time per energy
        numOfRays = size(obj.stf(i).ray, 2);
        for l = 1:numOfRays
            obj.stf(i).ray(l).energy = unique(obj.stf(i).ray(l).energy);
        end
        obj.stf(i).numOfRays = numel(obj.stf(i).ray);

        % save total number of bixels
        numOfBixels = 0;
        for j = 1:numel(obj.stf(i).ray)
            numOfBixels = numOfBixels + numel(obj.stf(i).ray(j).energy);
            obj.stf(i).numOfBixelsPerRay(j) = numel(obj.stf(i).ray(j).energy);
        end
        obj.stf(i).totalNumOfBixels = numOfBixels;

        % get bixelwidth
        bixelWidth_help = zeros(size(obj.stf(i).ray, 2), 2);
        for j = 1:obj.stf(i).numOfRays
            bixelWidth_help(j, 1) = obj.stf(i).ray(j).rayPos_bev(1);
            bixelWidth_help(j, 2) = obj.stf(i).ray(j).rayPos_bev(3);
        end
        bixelWidth_help1 = unique(round(1e3 * bixelWidth_help(:, 1)) / 1e3, 'sorted');
        bixelWidth_help2 = unique(round(1e3 * bixelWidth_help(:, 2)) / 1e3, 'sorted');

        bixelWidth = unique([unique(diff(bixelWidth_help1))' unique(diff(bixelWidth_help2))']);

        if numel(bixelWidth) == 1
            obj.stf(i).bixelWidth = bixelWidth;
        else
            obj.stf(i).bixelWidth = NaN;
        end
    end % deliveryType / ScanMode branch

end

%% check if matching given machine
% if a machine is given, check if we can exactly match the plan with
% the machine details
machineNotMatching = false;

if ~isempty(obj.pln.machine)
    matRad_cfg.dispInfo('Machine provided! Checking for exact steering match within RTPlan...');
    for i = 1:numel(obj.stf)
        for j = 1:obj.stf(i).numOfRays
            % loop over all energies
            numOfEnergy = length(obj.stf(i).ray(j).energy);
            for k = 1:numOfEnergy
                energyTemp = obj.stf(i).ray(j).energy(k);
                focusFWHM = obj.stf(i).ray(j).focusFWHM(k);
                energyIndex = find(abs([machine.data(:).energy] - energyTemp) < 10^-2);
                if isempty(energyIndex)
                    machineNotMatching = true;
                    break
                end
                focusIndex = find(abs([machine.data(energyIndex).initFocus.SisFWHMAtIso] - focusFWHM) < 10^-3);
                if isempty(focusIndex)
                    machineNotMatching = true;
                    break
                end
            end

            if machineNotMatching
                break
            end
        end
        if machineNotMatching
            break
        end
    end

    % If the machine matches, format the stf for direct use. Otherwise,
    % leave it be
    if machineNotMatching
        matRad_cfg.dispInfo('not matching!\n');
        matRad_cfg.dispWarning('The given machine does not match the steering info found in RTPlan. matRad will generate an stf, but it will be incompatible with the given machine and most likely not directly be usable in dose calculation!');

        for i = 1:numel(obj.stf)
            obj.stf(i).SAD = mean([obj.stf(i).SAD_x obj.stf(i).SAD_y]);
        end
    else
        matRad_cfg.dispInfo('matching!\n');
        matRad_cfg.dispInfo('Formatting stf for use with given machine...');

        for i = 1:numel(obj.stf)
            obj.stf(i).SAD = machine.meta.SAD;
            for j = 1:obj.stf(i).numOfRays

                % loop over all energies
                numOfEnergy = length(obj.stf(i).ray(j).energy);
                for k = 1:numOfEnergy
                    % If a corresponding machine was found, check assignment here
                    if ~isempty(machine)
                        energyTemp = obj.stf(i).ray(j).energy(k);
                        focusFWHM = obj.stf(i).ray(j).focusFWHM(k);
                        energyIndex = find(abs([machine.data(:).energy] - energyTemp) < 10^-2);
                        focusIndex = find(abs([machine.data(energyIndex).initFocus.SisFWHMAtIso] - focusFWHM) < 10^-3);

                        obj.stf(i).ray(j).energy(k) = machine.data(energyIndex).energy;
                        obj.stf(i).ray(j).focusIx(k) = focusIndex;
                        obj.stf(i).ray(j).focusFWHM(k) = machine.data(energyIndex).initFocus.SisFWHMAtIso(obj.stf(i).ray(j).focusIx(k));
                    end
                end
            end

        end
    end
end

%% Finalize geometry
for i = 1:numel(obj.stf)
    % coordinate transformation with rotation matrix.
    % use transpose matrix because we are working with row vectors
    rotMat_vectors_T = transpose(matRad_getRotationMatrix(obj.stf(i).gantryAngle, obj.stf(i).couchAngle));

    % set source point using (average/machine) SAD
    obj.stf(i).sourcePoint_bev = [0 -obj.stf(i).SAD 0];

    % Rotated Source point (1st gantry, 2nd couch)
    obj.stf(i).sourcePoint = obj.stf(i).sourcePoint_bev * rotMat_vectors_T;

    % Save ray and target position in lps system.
    for j = 1:obj.stf(i).numOfRays
        obj.stf(i).ray(j).targetPoint_bev = [2 * obj.stf(i).ray(j).rayPos_bev(1) obj.stf(i).SAD 2 * obj.stf(i).ray(j).rayPos_bev(3)];
        obj.stf(i).ray(j).rayPos      = obj.stf(i).ray(j).rayPos_bev * rotMat_vectors_T;
        obj.stf(i).ray(j).targetPoint = obj.stf(i).ray(j).targetPoint_bev * rotMat_vectors_T;
    end

    % book keeping & calculate focus index
    for j = 1:obj.stf(i).numOfRays
        obj.stf(i).numOfBixelsPerRay(j) = numel([obj.stf(i).ray(j).energy]);
    end

    obj.stf(i).timeStamp = datetime('now');
end

if any(isnan([obj.stf(:).bixelWidth])) || numel(unique([obj.stf(:).bixelWidth])) > 1
    obj.pln.propStf.bixelWidth = NaN;
else
    obj.pln.propStf.bixelWidth = obj.stf(1).bixelWidth;
end

end
