function matRad_makeTallyMCNP4biolOptiRMF(ct, pln, fileID_C_rest, binIntervals)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description goes here.
%
% call
%   matRad_makeTallyMCNP4biolOptiRMF(ct, pln, fileID_C_rest, binIntervals)
%
% input
%   ct, pln, fileID_C_rest, binIntervals
% output:
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User’s Manual. LACP-00634, May, 2013.
%   [2] Stewart et al., Phys. Med. Biol.,2015
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2020
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write/read parameters for repeated lattice tally
pln.propMCNP.latticeTally.cellLocal = ['((',int2str(2:size(binIntervals,2)), ')<', ...
    int2str(size(binIntervals,2)+1),'[0:', num2str(ct.cubeDim(2)-1), ' 0:', num2str(ct.cubeDim(1)-1), ' 0:', num2str(ct.cubeDim(3)-1), ']<',...
    int2str(size(binIntervals,2)+2), ')'];
pln.propMCNP.sdCardInfo = 1;

%% A. Read RBE DSB values
%% A.1 Identify source files
% Read RBE values for aerobic environment
fileNameRBEvalues.proton_aero = 'MCDS_protonRBE_DSB_aerobic*';
fileNameRBEvalues.deuteron_aero = 'MCDS_deuteronRBE_DSB_aerobic*';
fileNameRBEvalues.triton_aero = 'MCDS_tritonRBE_DSB_aerobic*';
fileNameRBEvalues.he3_aero = 'MCDS_3HeRBE_DSB_aerobic*';
fileNameRBEvalues.alpha_aero = 'MCDS_alphaRBE_DSB_aerobic*';
fileNameRBEvalues.electron_aero = 'MCDS_electronRBE_DSB_aerobic*';
fileNameRBEvalues.lithium_aero = 'MCDS_lithiumRBE_DSB_aerobic*';
% Read intra track values for aerobic environment
fileNameIntraTrackvalues.proton_aero = 'MCDS_protonIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.deuteron_aero = 'MCDS_deuteronIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.triton_aero = 'MCDS_tritonIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.he3_aero = 'MCDS_3HeIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.alpha_aero = 'MCDS_alphaIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.electron_aero = 'MCDS_electronIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.lithium_aero = 'MCDS_lithiumIntraTrack_intraT_aerobic*';
fileNameIntraTrackvalues.heavyIons_aero = 'MCDS_heavyIonsIntraTrack_intraT_aerobic*';

% Read RBE values for anoxic envirnmoent
fileNameRBEvalues.proton_anox = 'MCDS_protonRBE_DSB_anoxic*';
fileNameRBEvalues.deuteron_anox = 'MCDS_deuteronRBE_DSB_anoxic*';
fileNameRBEvalues.triton_anox = 'MCDS_tritonRBE_DSB_anoxic*';
fileNameRBEvalues.he3_anox = 'MCDS_3HeRBE_DSB_anoxic*';
fileNameRBEvalues.alpha_anox = 'MCDS_alphaRBE_DSB_anoxic*';
fileNameRBEvalues.electron_anox = 'MCDS_electronRBE_DSB_anoxic*';
fileNameRBEvalues.lithium_anox = 'MCDS_lithiumRBE_DSB_anoxic*';
% Read intra track values for anoxic envirnmoent
fileNameIntraTrackvalues.proton_anox = 'MCDS_protonIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.deuteron_anox = 'MCDS_deuteronIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.triton_anox = 'MCDS_tritonIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.he3_anox = 'MCDS_3HeIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.alpha_anox = 'MCDS_alphaIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.electron_anox = 'MCDS_electronIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.lithium_anox = 'MCDS_lithiumIntraTrack_intraT_anoxic*';
fileNameIntraTrackvalues.heavyIons_anox = 'MCDS_heavyIonsIntraTrack_intraT_anoxic*';

%% A.2 RBE values aerobic envirnment
% Read data from file: RBE factors for DSB in aerobic envirnmoent
RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
% Proton
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.proton_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
proton.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
RBEValues = rmfield(RBEValues, 'RBEValuesList');
% Deuteron
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.deuteron_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
deuteron.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
% Triton
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.triton_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
triton.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
% He3
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.he3_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
he3.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
% Alpha
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.alpha_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
alpha.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
% Electron
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.electron_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
electron.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);
% Lithium
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.lithium_aero]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
lithium.RBEValue_aero = RBEValue';
fclose(fid_RBEVal);

% Read data from file: Intra track term for RMF in aerobic envirnmoent
intraTrackValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
% Proton
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.proton_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
proton.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
intraTrackValues = rmfield(intraTrackValues, 'intraTrackValuesList');
% Deuteron
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.deuteron_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
deuteron.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% Triton
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.triton_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
triton.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% He3
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.he3_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
he3.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% Alpha
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.alpha_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
alpha.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% Electron
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.electron_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
electron.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% Lithium
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.lithium_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
lithium.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);
% Heavy ion
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.heavyIons_aero]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
heavyIons.intraTrackValue_aero = intraTrackValue';
fclose(fid_intraTrackVal);

%% A.3 RBE values anoxic envirnment
% Read data from file: RBE factors for DSB in anoxic envirnmoent
RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
% Proton
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.proton_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
proton.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
RBEValues = rmfield(RBEValues, 'RBEValuesList');
% Deuteron
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.deuteron_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
deuteron.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
% Triton
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.triton_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
triton.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
% He3
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.he3_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
he3.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
% Alpha
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.alpha_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
alpha.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
% Electron
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.electron_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
electron.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);
% Lithium
RBEValues.RBEValuesList = dir([RBEValues.pathLocation, fileNameRBEvalues.lithium_anox]);
fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
lithium.RBEValue_anox = RBEValue';
fclose(fid_RBEVal);

% Read data from file: Intra track term for RMF in anoxic envirnmoent
intraTrackValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
% Proton
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.proton_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
proton.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
intraTrackValues = rmfield(intraTrackValues, 'intraTrackValuesList');
% Deuteron
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.deuteron_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
deuteron.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% Triton
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.triton_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
triton.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% He3
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.he3_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
he3.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% Alpha
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.alpha_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
alpha.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% Electron
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.electron_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
electron.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% Lithium
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.lithium_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
lithium.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);
% Heavy ion
intraTrackValues.intraTrackValuesList = dir([intraTrackValues.pathLocation, fileNameIntraTrackvalues.heavyIons_anox]);
fid_intraTrackVal = fopen([intraTrackValues.intraTrackValuesList.folder, filesep, intraTrackValues.intraTrackValuesList.name], 'r');
intraTrackValue = fscanf(fid_intraTrackVal, '%f', [2,inf]);
heavyIons.intraTrackValue_anox = intraTrackValue';
fclose(fid_intraTrackVal);

%% B. Set up tallies: Dose and RBE calculation
%% B.1 Proton dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for proton dose calculation...')
disp('*****')
latticeTally.geometry = ['F1016:h ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Proton dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization (RBE*protonDose for aerobic environment)...')
disp('*****')
latticeTally.geometry = ['F1026:h ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(proton) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE1026\n');
for counterRBEvalues = 1:size(proton.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF1026\n');
for counterRBEvalues = 1:size(proton.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model
disp('*****')
disp('Tally type: Lattice tally for biological optimization (intra track term for aerobic environment)...')
disp('*****')
latticeTally.geometry = ['F1036:h ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(proton) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE1036\n');
for counterIntraTrackValues = 1:size(proton.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF1036\n');
for counterIntraTrackValues = 1:size(proton.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization (RBE*protonDose for anoxic environment)...')
disp('*****')
latticeTally.geometry = ['F1046:h ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(proton) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE1046\n');
for counterRBEvalues = 1:size(proton.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF1046\n');
for counterRBEvalues = 1:size(proton.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model
disp('*****')
disp('Tally type: Lattice tally for biological optimization (intra track term for anoxic environment)...')
disp('*****')
latticeTally.geometry = ['F1056:h ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(proton) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE1056\n');
for counterIntraTrackValues = 1:size(proton.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF1056\n');
for counterIntraTrackValues = 1:size(proton.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(proton.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.2 Deuteron dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for deuteron dose calculation...')
disp('*****')
latticeTally.geometry = ['F2016:d ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD2016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Deuteron dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*deuteronDose)...')
disp('*****')
latticeTally.geometry = ['F2026:d ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD2026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(deuteron) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE2026\n');
for counterRBEvalues = 1:size(deuteron.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF2026\n');
for counterRBEvalues = 1:size(deuteron.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F2036:d ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD2036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(deuteron) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE2036\n');
for counterIntraTrackValues = 1:size(deuteron.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF2036\n');
for counterIntraTrackValues = 1:size(deuteron.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization  for anoxic environment (RBE*deuteronDose)...')
disp('*****')
latticeTally.geometry = ['F2046:d ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD2046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(deuteron) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE2046\n');
for counterRBEvalues = 1:size(deuteron.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF2046\n');
for counterRBEvalues = 1:size(deuteron.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model
disp('*****')
disp('Tally type: Lattice tally for biological optimization  for anoxic environment  (intra track term)...')
disp('*****')
latticeTally.geometry = ['F2056:d ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD2056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(deuteron) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE2056\n');
for counterIntraTrackValues = 1:size(deuteron.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF2056\n');
for counterIntraTrackValues = 1:size(deuteron.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(deuteron.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.3 Triton dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for triton dose calculation...')
disp('*****')
latticeTally.geometry = ['F3016:t ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD3016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Triton dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*tritonDose)...')
disp('*****')
latticeTally.geometry = ['F3026:t ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD3026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(triton) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE3026\n');
for counterRBEvalues = 1:size(triton.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF3026\n');
for counterRBEvalues = 1:size(triton.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F3036:t ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD3036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(triton) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write triton DE/DF cards
fprintf(fileID_C_rest, 'DE3036\n');
for counterIntraTrackValues = 1:size(triton.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF3036\n');
for counterIntraTrackValues = 1:size(triton.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (RBE*tritonDose)...')
disp('*****')
latticeTally.geometry = ['F3046:t ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD3046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(triton) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE3046\n');
for counterRBEvalues = 1:size(triton.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF3046\n');
for counterRBEvalues = 1:size(triton.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F3056:t ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD3056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(triton) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write triton DE/DF cards
fprintf(fileID_C_rest, 'DE3056\n');
for counterIntraTrackValues = 1:size(triton.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF3056\n');
for counterIntraTrackValues = 1:size(triton.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(triton.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.4 He3 dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for He3 dose calculation...')
disp('*****')
latticeTally.geometry = ['F4016:s ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD4016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: He3Dose dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*He3Dose)...')
disp('*****')
latticeTally.geometry = ['F4026:s ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD4026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(He3Dose) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE4026\n');
for counterRBEvalues = 1:size(he3.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF4026\n');
for counterRBEvalues = 1:size(he3.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F4036:s ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD4036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(he3) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write He3 DE/DF cards
fprintf(fileID_C_rest, 'DE4036\n');
for counterIntraTrackValues = 1:size(he3.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF4036\n');
for counterIntraTrackValues = 1:size(he3.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (RBE*He3Dose)...')
disp('*****')
latticeTally.geometry = ['F4046:s ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD4046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(He3Dose) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE4046\n');
for counterRBEvalues = 1:size(he3.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF4046\n');
for counterRBEvalues = 1:size(he3.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F4056:s ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD4056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(he3) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write He3 DE/DF cards
fprintf(fileID_C_rest, 'DE4056\n');
for counterIntraTrackValues = 1:size(he3.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF4056\n');
for counterIntraTrackValues = 1:size(he3.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(he3.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.5 Alpha dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for alpha dose calculation...')
disp('*****')
latticeTally.geometry = ['F5016:a ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD5016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Alpha dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*alphaDose)...')
disp('*****')
latticeTally.geometry = ['F5026:a ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD5026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(alphaDose) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE5026\n');
for counterRBEvalues = 1:size(alpha.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF5026\n');
for counterRBEvalues = 1:size(alpha.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F5036:a ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD5036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(alpha) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write alpha DE/DF cards
fprintf(fileID_C_rest, 'DE5036\n');
for counterIntraTrackValues = 1:size(alpha.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF5036\n');
for counterIntraTrackValues = 1:size(alpha.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (RBE*alphaDose)...')
disp('*****')
latticeTally.geometry = ['F5046:a ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD5046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(alphaDose) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write proton DE/DF cards
fprintf(fileID_C_rest, 'DE5046\n');
for counterRBEvalues = 1:size(alpha.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF5046\n');
for counterRBEvalues = 1:size(alpha.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F5056:a ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD5056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(alpha) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write alpha DE/DF cards
fprintf(fileID_C_rest, 'DE5056\n');
for counterIntraTrackValues = 1:size(alpha.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF5056\n');
for counterIntraTrackValues = 1:size(alpha.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(alpha.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end
%% B.6 Electron dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for electron dose calculation...')
disp('*****')
latticeTally.geometry = ['F6016:e ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD6016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Electron dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*electronDose)...')
disp('*****')
latticeTally.geometry = ['F6026:e ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD6026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(electronDose) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE6026\n');
for counterRBEvalues = 1:size(electron.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF6026\n');
for counterRBEvalues = 1:size(electron.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F6036:e ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD6036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(electron) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE6036\n');
for counterIntraTrackValues = 1:size(electron.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF6036\n');
for counterIntraTrackValues = 1:size(electron.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (RBE*electronDose)...')
disp('*****')
latticeTally.geometry = ['F6046:e ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD6046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(electronDose) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE6046\n');
for counterRBEvalues = 1:size(electron.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF6046\n');
for counterRBEvalues = 1:size(electron.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F6056:e ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD6056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(electron) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE6056\n');
for counterIntraTrackValues = 1:size(electron.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF6056\n');
for counterIntraTrackValues = 1:size(electron.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(electron.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.7 Li-7 dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for lithium dose calculation...')
disp('*****')
latticeTally.geometry = ['F7016:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD7016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];
latticeTally.particleExtractionCard = 'FT7016 res 3007\n';

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Lithium dose for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
fprintf(fileID_C_rest, latticeTally.particleExtractionCard);

% RBExDose for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (RBE*lithiumDose)...')
disp('*****')
latticeTally.geometry = ['F7026:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD7026 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];
latticeTally.particleExtractionCard = 'FT7026 res 3007\n';

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(lithiumDose) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
fprintf(fileID_C_rest, latticeTally.particleExtractionCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE7026\n');
for counterRBEvalues = 1:size(lithium.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.RBEValue_aero(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF7026\n');
for counterRBEvalues = 1:size(lithium.RBEValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.RBEValue_aero(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F7036:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD7036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];
latticeTally.particleExtractionCard = 'FT7036 res 3007\n';

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(lithium) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
fprintf(fileID_C_rest, latticeTally.particleExtractionCard);

% Write electron DE/DF cards
fprintf(fileID_C_rest, 'DE7036\n');
for counterIntraTrackValues = 1:size(lithium.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF7036\n');
for counterIntraTrackValues = 1:size(lithium.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% RBExDose for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (RBE*lithiumDose)...')
disp('*****')
latticeTally.geometry = ['F7046:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD7046 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];
latticeTally.particleExtractionCard = 'FT7046 res 3007\n';

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: RBExDose(lithiumDose) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
fprintf(fileID_C_rest, latticeTally.particleExtractionCard);

% Write lithium DE/DF cards
fprintf(fileID_C_rest, 'DE7046\n');
for counterRBEvalues = 1:size(lithium.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.RBEValue_anox(counterRBEvalues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF7046\n');
for counterRBEvalues = 1:size(lithium.RBEValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.RBEValue_anox(counterRBEvalues,2)), '\n']);
end

% Intra Track Term RMF model for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F7056:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD7056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];
latticeTally.particleExtractionCard = 'FT7056 res 3007\n';

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(lithium) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
fprintf(fileID_C_rest, latticeTally.particleExtractionCard);

% Write lithium DE/DF cards
fprintf(fileID_C_rest, 'DE7056\n');
for counterIntraTrackValues = 1:size(lithium.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF7056\n');
for counterIntraTrackValues = 1:size(lithium.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(lithium.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.8 Heavy ion dose and RBE calculation
% Dose
disp('*****')
disp('Tally type: Lattice tally for heavy ion dose calculation...')
disp('*****')
latticeTally.geometry = ['F8016:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD8016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Heavy ion dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Intra Track Term RMF model for aerobic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for aerobic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F8036:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD8036 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(heavyIon) for aerobic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write heavy ion DE/DF cards
fprintf(fileID_C_rest, 'DE8036\n');
for counterIntraTrackValues = 1:size(heavyIons.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(heavyIons.intraTrackValue_aero(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF8036\n');
for counterIntraTrackValues = 1:size(heavyIons.intraTrackValue_aero,1)
    fprintf(fileID_C_rest, ['        ', num2str(heavyIons.intraTrackValue_aero(counterIntraTrackValues,2)), '\n']);
end

% Intra Track Term RMF model for anoxic environment
disp('*****')
disp('Tally type: Lattice tally for biological optimization for anoxic environment (intra track term)...')
disp('*****')
latticeTally.geometry = ['F8056:# ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD8056 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Intra Track Term(heavyIon) for anoxic environment (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

% Write heavy ion DE/DF cards
fprintf(fileID_C_rest, 'DE8056\n');
for counterIntraTrackValues = 1:size(heavyIons.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(heavyIons.intraTrackValue_anox(counterIntraTrackValues,1)), '\n']);
end
fprintf(fileID_C_rest, 'DF8056\n');
for counterIntraTrackValues = 1:size(heavyIons.intraTrackValue_anox,1)
    fprintf(fileID_C_rest, ['        ', num2str(heavyIons.intraTrackValue_anox(counterIntraTrackValues,2)), '\n']);
end

%% B.9 Neutron heating
disp('*****')
disp('Tally type: Lattice tally for neutron dose calculation...')
disp('*****')
latticeTally.geometry = ['F9016:n ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD9016 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Neutron dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);

%% B.10 Photon heating
disp('*****')
disp('Tally type: Lattice tally for photon dose calculation...')
disp('*****')
latticeTally.geometry = ['F1116:p ', pln.propMCNP.latticeTally.cellLocal, '\n'];
latticeTally.cellVolumeCard = ['SD1116 ', num2str(pln.propMCNP.sdCardInfo), ' ', num2str(pln.numOfVoxels-1), 'R\n'];

%Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Photon dose (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, latticeTally.geometry);
fprintf(fileID_C_rest, latticeTally.cellVolumeCard);
