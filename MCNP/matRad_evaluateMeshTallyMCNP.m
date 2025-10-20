function doseMatrix = matRad_evaluateMeshTallyMCNP(fileName, pln, ct)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate files containing tally results of MCNP calculation.
%
% call
%   doseMatrix = matRad_evaluateMeshTallyMCNP(fileName, tallyType, ctData)
%
% input
%   fileName:       file name of the text file to evaluate
%   pln.propMCNP.tallySpecifier:    type of tally used for simulation
%                                   possible input: 'KERMA_F4' or
%                                   'TotalDose_TMESH'
%   ctData:         3D volume from ct data containing tissue
%                   characteristics for every voxel
%
%
% output
%   doseMatrix:     3D matrix containing result of MCNP calculation
%                   .neutronDose: KERMA from neutron interaction
%                   (F4 tally)
%                   .photonDose: KERMA from photon interaction
%                   (F4 tally)
%                   .dose: dose as total depositied energy in material
%                   (TMESH tally)
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 10/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TMESH Tally for Total Dose
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TMESH tally detects all energy deposited in the medium by taking account
% also of secondary particles (if not switched off on mode card).
%% Read MCNP result for TMESH type 3 tally
%     fid_read = fopen(fileName, 'r');
%     dummyLine = 'dummyLine'; % Use a dummy line to find data values
%     while ~strcmp(dummyLine, 'vals'), dummyLine = fgetl(fid_read); end % Set location in mctal file to position where values start

%     resultMCNP_dummy = zeros(prod(ct.cubeDim)*2,1)';
%     resultMCNP_dummy = fscanf(fid_read, '%f', [1 inf]); % Values are given with relative errors
%
%     resultMCNP(:,1) = resultMCNP_dummy(1:2:end);
%     resultMCNP(:,2) = resultMCNP_dummy(2:2:end);

resultMCNP = matRad_readDataFromText_TMESHvBioOpti(fileName, 'TMESH3', 2);
resultMCNP = resultMCNP';
doseMatrix.physicalDose = zeros(ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3));   % Total dose
doseMatrix.physicalDose(1:end) = resultMCNP(:,1);


doseMatrix.physicalDose_relError = zeros(ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3));  % Relative error of total dose
doseMatrix.physicalDose_relError(1:end) = resultMCNP(:,2);

clear resultMCNP
%     fclose(fid_read);

doseMatrix.physicalDose = permute(doseMatrix.physicalDose, [2,1,3]);    % Permute matrix to match matRad coordinate system
for cellCounter = 1:size(ct.tissueBin,2)
    if cellCounter == 1
        doseMatrix.physicalDose(ct.tissueBin(cellCounter).linIndVol) = 0;   % Dose deposition in air is neglected
    else
        doseMatrix.physicalDose(ct.tissueBin(cellCounter).linIndVol) = ...
            doseMatrix.physicalDose(ct.tissueBin(cellCounter).linIndVol)./mean(ct.density{1,1}(ct.tissueBin(cellCounter).linIndVol));   % MCNP output is in MeV/cm^3/source particle & ct.density is given in g/cm^3
    end
    mean(ct.density{1,1}(ct.tissueBin(cellCounter).linIndVol))
end
doseMatrix.physicalDose = doseMatrix.physicalDose*1.602177e-19*1e6*1e3; % Convert MeV/g to J/kg, output is now in Gy/source particle

doseMatrix.physicalDose_relError = permute(doseMatrix.physicalDose_relError, [2,1,3]);

% Bad TMESH tally statistics can lead to negative results
if sum(doseMatrix.physicalDose<0, 'all')
    warning('*********************')
    warning('Negative TMESH tally results detected. This is a hint for bad statistics!')
    warning(['Minimum value: ', num2str(min(doseMatrix.physicalDose, [], 'all')), ' Gy and maximum value: ', num2str(max(doseMatrix.physicalDose, [], 'all')), ' Gy.'])
    warning('*********************')

    doseMatrix.physicalDose(doseMatrix.physicalDose<0) = 0;
    doseMatrix.physicalDose_relError(doseMatrix.physicalDose<0) = 0;
end
%% Read RBE weighted dose
if isfield(pln.propOpt,'bioOptimization')
    if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
        % Read data
        tallyData = matRad_readDataFromText_TMESHvBioOpti(fileName, 'F6heating4RBEcalc', 2);

        % Read proton data
        %              [doseMatrix.physicalDose_proton,  doseMatrix.physicalDose_proton_relError] = convertDose4RBECalculation(ct, tallyData.tally1016);
        %              [doseMatrix.RBExDose_proton_aero, doseMatrix.RBExDose_proton_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally1026);
        %              [doseMatrix.intraTrackTerm_aero_proton, doseMatrix.intraTrackTerm_proton_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally1036);
        %              [doseMatrix.RBExDose_proton_anox, doseMatrix.RBExDose_proton_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally1046);
        %              [doseMatrix.intraTrackTerm_anox_proton, doseMatrix.intraTrackTerm_proton_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally1056);

        % Read deuteron data
        %              [doseMatrix.physicalDose_deuteron,  doseMatrix.physicalDose_deuteron_relError] = convertDose4RBECalculation(ct, tallyData.tally2016);
        %              [doseMatrix.RBExDose_deuteron_aero, doseMatrix.RBExDose_deuteron_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally2026);
        %              [doseMatrix.intraTrackTerm_aero_deuteron, doseMatrix.intraTrackTerm_deuteron_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally2036);
        %             [doseMatrix.RBExDose_deuteron_anox, doseMatrix.RBExDose_deuteron_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally2046);
        %             [doseMatrix.intraTrackTerm_anox_deuteron, doseMatrix.intraTrackTerm_deuteron_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally2056);

        % Read triton data
        %              [doseMatrix.physicalDose_triton,  doseMatrix.physicalDose_triton_relError] = convertDose4RBECalculation(ct, tallyData.tally3016);
        %              [doseMatrix.RBExDose_triton_aero, doseMatrix.RBExDose_triton_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally3026);
        %              [doseMatrix.intraTrackTerm_aero_triton, doseMatrix.intraTrackTerm_triton_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally3036);
        %             [doseMatrix.RBExDose_triton_anox, doseMatrix.RBExDose_triton_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally3046);
        %             [doseMatrix.intraTrackTerm_anox_triton, doseMatrix.intraTrackTerm_triton_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally3056);

        % Read He3 data
        %             [doseMatrix.physicalDose_he3,  doseMatrix.physicalDose_he3_relError] = convertDose4RBECalculation(ct, tallyData.tally4016);
        %             [doseMatrix.RBExDose_he3_aero, doseMatrix.RBExDose_he3_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally4026);
        %             [doseMatrix.intraTrackTerm_aero_he3, doseMatrix.intraTrackTerm_he3_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally4036);
        %             [doseMatrix.RBExDose_he3_anox, doseMatrix.RBExDose_he3_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally4046);
        %             [doseMatrix.intraTrackTerm_anox_he3, doseMatrix.intraTrackTerm_he3_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally4056);

        % Read alpha data
        %             [doseMatrix.physicalDose_alpha,  doseMatrix.physicalDose_alpha_relError] = convertDose4RBECalculation(ct, tallyData.tally5016);
        %             [doseMatrix.RBExDose_alpha_aero, doseMatrix.RBExDose_alpha_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally5026);
        %             [doseMatrix.intraTrackTerm_aero_alpha, doseMatrix.intraTrackTerm_alpha_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally5036);
        %             [doseMatrix.RBExDose_alpha_anox, doseMatrix.RBExDose_alpha_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally5046);
        %             [doseMatrix.intraTrackTerm_anox_alpha, doseMatrix.intraTrackTerm_alpha_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally5056);

        % Read electron data
        %             [doseMatrix.physicalDose_electron,  doseMatrix.physicalDose_electron_relError] = convertDose4RBECalculation(ct, tallyData.tally6016);
        %             [doseMatrix.RBExDose_electron_aero, doseMatrix.RBExDose_electron_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally6026);
        %             [doseMatrix.intraTrackTerm_aero_electron, doseMatrix.intraTrackTerm_electron_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally6036);
        %             [doseMatrix.RBExDose_electron_anox, doseMatrix.RBExDose_electron_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally6046);
        %             [doseMatrix.intraTrackTerm_anox_electron, doseMatrix.intraTrackTerm_electron_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally6056);

        % Read lithium data
        %             [doseMatrix.physicalDose_lithium,  doseMatrix.physicalDose_lithium_relError] = convertDose4RBECalculation(ct, tallyData.tally7016(:,1:size(tallyData.tally7016,2)/2));
        %             [doseMatrix.RBExDose_lithium_aero, doseMatrix.RBExDose_lithium_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally7026(:,1:size(tallyData.tally7026,2)/2));
        %             [doseMatrix.intraTrackTerm_aero_lithium, doseMatrix.intraTrackTerm_lithium_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally7036(:,1:size(tallyData.tally7036,2)/2));
        %             [doseMatrix.RBExDose_lithium_anox, doseMatrix.RBExDose_lithium_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally7046(:,1:size(tallyData.tally7046,2)/2));
        %             [doseMatrix.intraTrackTerm_anox_lithium, doseMatrix.intraTrackTerm_lithium_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally7056(:,1:size(tallyData.tally7056,2)/2));

        % Read heavy ion data
        %             [doseMatrix.physicalDose_heavyIon,  doseMatrix.physicalDose_heavyIon_relError] = convertDose4RBECalculation(ct, tallyData.tally8016);
        %             [doseMatrix.intraTrackTerm_aero_heavyIon, doseMatrix.intraTrackTerm_heavyIon_aero_relError] = convertDose4RBECalculation(ct, tallyData.tally8036);
        %             [doseMatrix.intraTrackTerm_anox_heavyIon, doseMatrix.intraTrackTerm_heavyIon_anox_relError] = convertDose4RBECalculation(ct, tallyData.tally8056);


        % Read neutron heating data
        %             [doseMatrix.physicalDose_neutronHeating, doseMatrix.physicalDose_neutronHeating_relError] = convertDose4RBECalculation(ct, tallyData.tally9016);

        % Read photon heating data
        %             [doseMatrix.physicalDose_photonHeating, doseMatrix.physicalDose_photonHeating_relError] = convertDose4RBECalculation(ct, tallyData.tally1116);
    end
end
%% Function Definition: matRad_cart2linIndex gives back linear indices of tally cells/voxels
    function linearIndex_Tally = matRad_cart2linIndex(resultMCNP,ct)

        grid_matInd = zeros(prod(ct.cubeDim),3);

        % MATLAB matrix indexing has to be considered: x-coordinate = column index
        grid_matInd(:,1) = round((resultMCNP(:,2)+(ct.resolution.y/2))/ct.resolution.y); % row
        grid_matInd(:,2) = round((resultMCNP(:,1)+(ct.resolution.x/2))/ct.resolution.x); % column
        grid_matInd(:,3) = round((resultMCNP(:,3)+(ct.resolution.z/2))/ct.resolution.z);

        % Calculate linear indices
        linearIndex_Tally = sub2ind(ct.cubeDim, grid_matInd(:,1), grid_matInd(:,2), grid_matInd(:,3));
    end

%% Function definition to read physical and RBE weighted dose for secondary particles
    function [physicalDose, physicalDose_relError] = convertDose4RBECalculation(ct, resultData)
        resultData = resultData';

        physicalDose = zeros(ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3));   % Total dose
        physicalDose(1:end) = resultData(:,1);


        physicalDose_relError = zeros(ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3));  % Relative error of total dose
        physicalDose_relError(1:end) = resultData(:,2);

        physicalDose = permute(physicalDose, [2,1,3]);    % Permute matrix to match matRad coordinate system

        % Conversion of tally data [MeV/g] where cell mass was set to 1g to Gy/source particle
        dummyVoxelVolume = ct.resolution.x*ct.resolution.y*ct.resolution.z*1e-3; % voxel volume in ccm

        for cellCounter_subFunc = 1:size(ct.tissueBin,2)
            if cellCounter_subFunc == 1 || cellCounter_subFunc == 2
                physicalDose(ct.tissueBin(cellCounter_subFunc).linIndVol) = 0;   % Dose deposition in air is neglected
            else
                physicalDose(ct.tissueBin(cellCounter_subFunc).linIndVol) = ...
                    physicalDose(ct.tissueBin(cellCounter_subFunc).linIndVol)./(mean(ct.density{1,1}(ct.tissueBin(cellCounter_subFunc).linIndVol))*dummyVoxelVolume);   % MCNP output is in MeV/cm^3/source particle & ct.density is given in g/cm^3
            end
            mean(ct.density{1,1}(ct.tissueBin(cellCounter_subFunc).linIndVol));
        end

        physicalDose = physicalDose*1.602177e-19*1e6*1e3; % Convert MeV/g to J/kg, output is now in Gy/source particle

        physicalDose_relError = permute(physicalDose_relError, [2,1,3]);
    end
end