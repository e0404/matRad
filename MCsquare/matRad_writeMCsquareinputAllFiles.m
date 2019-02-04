function matRad_writeMCsquareinputAllFiles(filename,MCsquareConfig,stf)
% generate input files for MCsquare dose calcualtion from matRad
% 
% call
%   matRad_writeMCsquareinputAllFiles(filename,MCsquareConfig,stf)
%
% input
%   filename:       filename
%   MCsquareConfig: matRad MCsquare configuration
%   stf:            matRad steering information struct
%
% output
%   -
%
% References
%   [1] https://openreggui.org/git/open/REGGUI/blob/master/functions/io/convert_Plan_PBS.m
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% write overall configuration file
fileHandle = fopen(filename,'w');
MCsquareConfig.write(fileHandle);
fclose(fileHandle);

%% prepare steering file writing
numOfFields = length(stf);
if MCsquareConfig.Beamlet_Mode
    totalMetersetWeightOfAllFields = 1;
else
     totalMetersetWeightOfFields = NaN*ones(numOfFields,1);
     for i = 1:numOfFields
        totalMetersetWeightOfFields(i) = sum([stf(i).energyLayer.numOfPrimaries]);
    end
    totalMetersetWeightOfAllFields = sum(totalMetersetWeightOfFields);
end

%% write steering file

fileHandle = fopen(MCsquareConfig.BDL_Plan_File,'w');

fprintf(fileHandle,'#TREATMENT-PLAN-DESCRIPTION\n');
fprintf(fileHandle,'#PlanName\n');
fprintf(fileHandle,'matRad_bixel\n');
fprintf(fileHandle,'#NumberOfFractions\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'##FractionID\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'##NumberOfFields\n');
fprintf(fileHandle,[num2str(numOfFields) '\n']);
for i = 1:numOfFields
    fprintf(fileHandle,'###FieldsID\n');
    fprintf(fileHandle,[num2str(i) '\n']);
end
fprintf(fileHandle,'\n#TotalMetersetWeightOfAllFields\n');
fprintf(fileHandle,[num2str(totalMetersetWeightOfAllFields) '\n']);
    
for i = 1:numOfFields
    fprintf(fileHandle,'\n#FIELD-DESCRIPTION\n');
    fprintf(fileHandle,'###FieldID\n');
    fprintf(fileHandle,[num2str(i) '\n']);
    fprintf(fileHandle,'###FinalCumulativeMeterSetWeight\n');
    if MCsquareConfig.Beamlet_Mode
        finalCumulativeMeterSetWeight = 1/numOfFields;
    else
        finalCumulativeMeterSetWeight = totalMetersetWeightOfFields(i);
    end
    fprintf(fileHandle,[num2str(finalCumulativeMeterSetWeight) '\n']);
    fprintf(fileHandle,'###GantryAngle\n');
    fprintf(fileHandle,[num2str(stf(i).gantryAngle) '\n']);
    fprintf(fileHandle,'###PatientSupportAngle\n');
    fprintf(fileHandle,[num2str(stf(i).couchAngle) '\n']);
    fprintf(fileHandle,'###IsocenterPosition\n');
    fprintf(fileHandle,[num2str(stf(i).isoCenter) '\n']);
    fprintf(fileHandle,'###NumberOfControlPoints\n');
    numOfEnergies = numel(stf(i).energies);
    fprintf(fileHandle,[num2str(numOfEnergies) '\n']);

    metersetOffset = 0;
    fprintf(fileHandle,'\n#SPOTS-DESCRIPTION\n');
    for j = 1:numOfEnergies
        fprintf(fileHandle,'####ControlPointIndex\n');
        fprintf(fileHandle,[num2str(j) '\n']);
        fprintf(fileHandle,'####SpotTunnedID\n');
        fprintf(fileHandle,['1\n']);
        fprintf(fileHandle,'####CumulativeMetersetWeight\n');
        if MCsquareConfig.Beamlet_Mode
            cumulativeMetersetWeight = j/numOfEnergies * 1/numOfFields;
        else 
            cumulativeMetersetWeight = metersetOffset + sum([stf(i).energyLayer(j).numOfPrimaries]);
            metersetOffset = cumulativeMetersetWeight;
        end
        fprintf(fileHandle,[num2str(cumulativeMetersetWeight) '\n']);
        fprintf(fileHandle,'####Energy (MeV)\n');
        fprintf(fileHandle,[num2str(stf(i).energies(j)) '\n']);
        fprintf(fileHandle,'####NbOfScannedSpots\n');
        numOfSpots = size(stf(i).energyLayer(j).targetPoints,1);
        fprintf(fileHandle,[num2str(numOfSpots) '\n']);
        fprintf(fileHandle,'####X Y Weight\n');
        for k = 1:numOfSpots
            if MCsquareConfig.Beamlet_Mode
                n = stf(i).energyLayer(j).numOfPrimaries(k);
            else
                n = stf(i).energyLayer(j).numOfPrimaries(k) / mcSquare_magicFudge(stf(i).energies(j));
            end
            fprintf(fileHandle,[num2str(stf(i).energyLayer(j).targetPoints(k,:)) ' ' num2str(n) '\n']);
        end
    end        
end

fclose(fileHandle);

end

function gain = mcSquare_magicFudge(energy)
% mcSquare will scale the spot intensities in
% https://gitlab.com/openmcsquare/MCsquare/blob/master/src/data_beam_model.c#L906
% by this factor so we need to divide up front to make things work. The
% original code can be found at https://gitlab.com/openmcsquare/MCsquare/blob/master/src/compute_beam_model.c#L16

K = 35.87; % in eV (other value 34.23 ?)

% // Air stopping power (fit ICRU) multiplied by air density
SP = (9.6139e-9*energy^4 - 7.0508e-6*energy^3 + 2.0028e-3*energy^2 - 2.7615e-1*energy + 2.0082e1) * 1.20479E-3 * 1E6; % // in eV / cm

% // Temp & Pressure correction
PTP = 1.0; 

% // MU calibration (1 MU = 3 nC/cm)
% // 1cm de gap effectif
C = 3.0E-9; % // in C / cm

% // Gain: 1eV = 1.602176E-19 J
gain = (C*K) / (SP*PTP*1.602176E-19);

% divide by 1e7 to not get tiny numbers...
gain = gain/1e7;

end