function matRad_writeMCsquareinputFiles(filename,MCsquareConfig,stf)


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
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'###FieldsID\n');
fprintf(fileHandle,'1\n\n');

fprintf(fileHandle,'#TotalMetersetWeightOfAllFields\n');
fprintf(fileHandle,'1000\n');
 
fprintf(fileHandle,'#FIELD-DESCRIPTION\n');
fprintf(fileHandle,'###FieldID\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'###FinalCumulativeMeterSetWeight\n');
fprintf(fileHandle,'1000\n');
fprintf(fileHandle,'###GantryAngle\n');
fprintf(fileHandle,[num2str(stf.gantryAngle) '\n']);
fprintf(fileHandle,'###PatientSupportAngle\n');
fprintf(fileHandle,[num2str(stf.couchAngle) '\n']);
fprintf(fileHandle,'###IsocenterPosition\n');
fprintf(fileHandle,[num2str(stf.isoCenter) '\n']);
fprintf(fileHandle,'###NumberOfControlPoints\n');
fprintf(fileHandle,'1\n');

fprintf(fileHandle,'#SPOTS-DESCRIPTION\n');
fprintf(fileHandle,'####ControlPointIndex\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'####SpotTunnedID\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'####CumulativeMetersetWeight\n');
fprintf(fileHandle,'1000\n');
fprintf(fileHandle,'####Energy (MeV)\n');
fprintf(fileHandle,[num2str(stf.energy) '\n']);
fprintf(fileHandle,'####NbOfScannedSpots\n');
fprintf(fileHandle,'1\n');
fprintf(fileHandle,'####X Y Weight\n');
fprintf(fileHandle,[num2str(stf.targetPoint) ' ' num2str(MCsquareConfig.Num_Primaries) '\n']);

fclose(fileHandle);


