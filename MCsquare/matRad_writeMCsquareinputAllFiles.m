function matRad_writeMCsquareinputAllFiles(filename,MCsquareConfig,stf)


%% write overall configuration file
fileHandle = fopen(filename,'w');
MCsquareConfig.write(fileHandle);
fclose(fileHandle);


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
numOfFields = length(stf);
fprintf(fileHandle,[num2str(numOfFields) '\n']);
for i = 1:numOfFields
    fprintf(fileHandle,'###FieldsID\n');
    fprintf(fileHandle,[num2str(i) '\n\n']);

    fprintf(fileHandle,'#TotalMetersetWeightOfAllFields\n');
    fprintf(fileHandle,'1000\n');

    fprintf(fileHandle,'#FIELD-DESCRIPTION\n');
    fprintf(fileHandle,'###FieldID\n');
    fprintf(fileHandle,[num2str(i) '\n']);
    fprintf(fileHandle,'###FinalCumulativeMeterSetWeight\n');
    fprintf(fileHandle,'1000\n');
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
        fprintf(fileHandle,'1000\n');
        fprintf(fileHandle,'####Energy (MeV)\n');
        fprintf(fileHandle,[num2str(stf(i).energies(j)) '\n']);
        fprintf(fileHandle,'####NbOfScannedSpots\n');
        numOfSpots = size(stf(i).energyLayer(j).targetPoints,1);
        fprintf(fileHandle,[num2str(numOfSpots) '\n']);
        fprintf(fileHandle,'####X Y Weight\n');
        for k = 1:numOfSpots
            fprintf(fileHandle,[num2str(stf(i).energyLayer(j).targetPoints(k,:)) ' ' num2str(stf(i).energyLayer(j).numOfPrimaries(k)) '\n']);
        end
    end        
end

fclose(fileHandle);


