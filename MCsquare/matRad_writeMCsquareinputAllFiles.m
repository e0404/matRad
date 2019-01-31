function matRad_writeMCsquareinputAllFiles(filename,MCsquareConfig,stf)


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
fprintf(fileHandle,[num2str(totalMetersetWeightOfAllFields) '\n\n']);
    
for i = 1:numOfFields
    fprintf(fileHandle,'#FIELD-DESCRIPTION\n');
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
    for j = 1:numOfEnergies
        fprintf(fileHandle,'#SPOTS-DESCRIPTION\n');
        fprintf(fileHandle,'####ControlPointIndex\n');
        fprintf(fileHandle,[num2str(j) '\n']);
        fprintf(fileHandle,'####SpotTunnedID\n');
        fprintf(fileHandle,[num2str(j) '\n']);
        fprintf(fileHandle,'####CumulativeMetersetWeight\n');
        if MCsquareConfig.Beamlet_Mode
            cumulativeMetersetWeight = 1/numOfEnergies * 1/numOfFields;
        else 
            cumulativeMetersetWeight = sum([stf(i).energyLayer(j).numOfPrimaries]);
        end
        fprintf(fileHandle,[num2str(cumulativeMetersetWeight) '\n']);
        fprintf(fileHandle,'####Energy (MeV)\n');
        fprintf(fileHandle,[num2str(stf(i).energies(j)) '\n']);
        fprintf(fileHandle,'####NbOfScannedSpots\n');
        numOfSpots = size(stf(i).energyLayer(j).targetPoints,1);
        fprintf(fileHandle,[num2str(numOfSpots) '\n']);
        fprintf(fileHandle,'####X Y Weight\n');
        for k = 1:numOfSpots
            n = stf(i).energyLayer(j).numOfPrimaries(k);
            fprintf(fileHandle,[num2str(stf(i).energyLayer(j).targetPoints(k,:)) ' ' num2str(n) '\n']);
        end
    end        
end

fclose(fileHandle);


