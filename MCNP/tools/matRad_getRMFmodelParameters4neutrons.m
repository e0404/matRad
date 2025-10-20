function [dij, cst] = matRad_getRMFmodelParameters4neutrons(dij, counterDijColumns, doseMatrixBixel, pln, ct, cst, modificationMode)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A: Extract data tallied in MCNP simulation when biological optimization
% for neutrons is switched on.
%
% B: Calculate RBE_LD and RBE_HD values for neutron radiation for each
% bixel. 
%
% input:        dij:                structure from matRad
%               counterDijColumns:  indicates bixel number
%               doseMatrixBixel:    information on data tallied in MCNP
%                                   calculation (only needed in simulation
%                                   workflow, set doseMatrixBixel=0 when 
%                                   not needed)
%               pln:                structure from matRad
%               ct:                 structure from matRad
%               cst:                structure from matRad (modify
%                                   alpha/beta ratio here)
%               modificationMode:   indicate modification of alpha/beta
%                                   ratio for RMF model parameter
%                                   re-calculation
%
%
% output:       dij
%               cst:                cst with cancelled alpha/beta ratio
%                                   overlaps
%
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 09/2020
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check function usage
if ~exist('modificationMode')
    modificationMode = false;
end

%% A: Data Extraction: Extract particle specific dose, RBE_DSB weighted dose and intra track term for RMF model
if ~modificationMode
    %% Proton
    % Physical dose
    dij.doseMatrixBixel(counterDijColumns).protonDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_proton),1,doseMatrixBixel.physicalDose_proton(doseMatrixBixel.physicalDose_proton~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).protonDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_proton_relError),1,doseMatrixBixel.physicalDose_proton_relError(doseMatrixBixel.physicalDose_proton_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).protonDosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_proton_aero),1,doseMatrixBixel.RBExDose_proton_aero(doseMatrixBixel.RBExDose_proton_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_proton_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_proton_aero_relError),1,doseMatrixBixel.RBExDose_proton_aero_relError(doseMatrixBixel.RBExDose_proton_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_proton{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_proton),1,doseMatrixBixel.intraTrackTerm_aero_proton(doseMatrixBixel.intraTrackTerm_aero_proton~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_proton_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_proton_aero_relError),1,doseMatrixBixel.intraTrackTerm_proton_aero_relError(doseMatrixBixel.intraTrackTerm_proton_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).protonDosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_proton_anox),1,doseMatrixBixel.RBExDose_proton_anox(doseMatrixBixel.RBExDose_proton_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_proton_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_proton_anox_relError),1,doseMatrixBixel.RBExDose_proton_anox_relError(doseMatrixBixel.RBExDose_proton_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_proton{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_proton),1,doseMatrixBixel.intraTrackTerm_anox_proton(doseMatrixBixel.intraTrackTerm_anox_proton~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_proton_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_proton_anox_relError),1,doseMatrixBixel.intraTrackTerm_proton_anox_relError(doseMatrixBixel.intraTrackTerm_proton_anox_relError~=0),dij.numOfVoxels,1);

    
    %% Deuteron
    dij.doseMatrixBixel(counterDijColumns).deuteronDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_deuteron),1,doseMatrixBixel.physicalDose_deuteron(doseMatrixBixel.physicalDose_deuteron~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).deuteronDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_deuteron_relError),1,doseMatrixBixel.physicalDose_deuteron_relError(doseMatrixBixel.physicalDose_deuteron_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).deuteronDosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_deuteron_aero),1,doseMatrixBixel.RBExDose_deuteron_aero(doseMatrixBixel.RBExDose_deuteron_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_deuteron_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_deuteron_aero_relError),1,doseMatrixBixel.RBExDose_deuteron_aero_relError(doseMatrixBixel.RBExDose_deuteron_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_deuteron{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_deuteron),1,doseMatrixBixel.intraTrackTerm_aero_deuteron(doseMatrixBixel.intraTrackTerm_aero_deuteron~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_deuteron_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_deuteron_aero_relError),1,doseMatrixBixel.intraTrackTerm_deuteron_aero_relError(doseMatrixBixel.intraTrackTerm_deuteron_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).deuteronDosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_deuteron_anox),1,doseMatrixBixel.RBExDose_deuteron_anox(doseMatrixBixel.RBExDose_deuteron_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_deuteron_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_deuteron_anox_relError),1,doseMatrixBixel.RBExDose_deuteron_anox_relError(doseMatrixBixel.RBExDose_deuteron_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_deuteron{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_deuteron),1,doseMatrixBixel.intraTrackTerm_anox_deuteron(doseMatrixBixel.intraTrackTerm_anox_deuteron~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_deuteron_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_deuteron_anox_relError),1,doseMatrixBixel.intraTrackTerm_deuteron_anox_relError(doseMatrixBixel.intraTrackTerm_deuteron_anox_relError~=0),dij.numOfVoxels,1);
    
    %% Triton
    dij.doseMatrixBixel(counterDijColumns).tritonDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_triton),1,doseMatrixBixel.physicalDose_triton(doseMatrixBixel.physicalDose_triton~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).tritonDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_triton_relError),1,doseMatrixBixel.physicalDose_triton_relError(doseMatrixBixel.physicalDose_triton_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).tritonDosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_triton_aero),1,doseMatrixBixel.RBExDose_triton_aero(doseMatrixBixel.RBExDose_triton_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_triton_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_triton_aero_relError),1,doseMatrixBixel.RBExDose_triton_aero_relError(doseMatrixBixel.RBExDose_triton_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_triton{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_triton),1,doseMatrixBixel.intraTrackTerm_aero_triton(doseMatrixBixel.intraTrackTerm_aero_triton~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_triton_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_triton_aero_relError),1,doseMatrixBixel.intraTrackTerm_triton_aero_relError(doseMatrixBixel.intraTrackTerm_triton_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).tritonDosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_triton_anox),1,doseMatrixBixel.RBExDose_triton_anox(doseMatrixBixel.RBExDose_triton_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_triton_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_triton_anox_relError),1,doseMatrixBixel.RBExDose_triton_anox_relError(doseMatrixBixel.RBExDose_triton_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_triton{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_triton),1,doseMatrixBixel.intraTrackTerm_anox_triton(doseMatrixBixel.intraTrackTerm_anox_triton~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_triton_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_triton_anox_relError),1,doseMatrixBixel.intraTrackTerm_triton_anox_relError(doseMatrixBixel.intraTrackTerm_triton_anox_relError~=0),dij.numOfVoxels,1);
        
    %% He3
    dij.doseMatrixBixel(counterDijColumns).he3Dose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_he3),1,doseMatrixBixel.physicalDose_he3(doseMatrixBixel.physicalDose_he3~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).he3Dose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_he3_relError),1,doseMatrixBixel.physicalDose_he3_relError(doseMatrixBixel.physicalDose_he3_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).he3DosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_he3_aero),1,doseMatrixBixel.RBExDose_he3_aero(doseMatrixBixel.RBExDose_he3_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_he3_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_he3_aero_relError),1,doseMatrixBixel.RBExDose_he3_aero_relError(doseMatrixBixel.RBExDose_he3_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_he3{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_he3),1,doseMatrixBixel.intraTrackTerm_aero_he3(doseMatrixBixel.intraTrackTerm_aero_he3~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_he3_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_he3_aero_relError),1,doseMatrixBixel.intraTrackTerm_he3_aero_relError(doseMatrixBixel.intraTrackTerm_he3_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).he3DosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_he3_anox),1,doseMatrixBixel.RBExDose_he3_anox(doseMatrixBixel.RBExDose_he3_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_he3_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_he3_anox_relError),1,doseMatrixBixel.RBExDose_he3_anox_relError(doseMatrixBixel.RBExDose_he3_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_he3{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_he3),1,doseMatrixBixel.intraTrackTerm_anox_he3(doseMatrixBixel.intraTrackTerm_anox_he3~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_he3_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_he3_anox_relError),1,doseMatrixBixel.intraTrackTerm_he3_anox_relError(doseMatrixBixel.intraTrackTerm_he3_anox_relError~=0),dij.numOfVoxels,1);
   
    %% Alpha
    dij.doseMatrixBixel(counterDijColumns).alphaDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_alpha),1,doseMatrixBixel.physicalDose_alpha(doseMatrixBixel.physicalDose_alpha~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).alphaDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_alpha_relError),1,doseMatrixBixel.physicalDose_alpha_relError(doseMatrixBixel.physicalDose_alpha_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).alphaDosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_alpha_aero),1,doseMatrixBixel.RBExDose_alpha_aero(doseMatrixBixel.RBExDose_alpha_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_alpha_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_alpha_aero_relError),1,doseMatrixBixel.RBExDose_alpha_aero_relError(doseMatrixBixel.RBExDose_alpha_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_alpha{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_alpha),1,doseMatrixBixel.intraTrackTerm_aero_alpha(doseMatrixBixel.intraTrackTerm_aero_alpha~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_alpha_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_alpha_aero_relError),1,doseMatrixBixel.intraTrackTerm_alpha_aero_relError(doseMatrixBixel.intraTrackTerm_alpha_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).alphaDosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_alpha_anox),1,doseMatrixBixel.RBExDose_alpha_anox(doseMatrixBixel.RBExDose_alpha_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_alpha_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_alpha_anox_relError),1,doseMatrixBixel.RBExDose_alpha_anox_relError(doseMatrixBixel.RBExDose_alpha_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_alpha{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_alpha),1,doseMatrixBixel.intraTrackTerm_anox_alpha(doseMatrixBixel.intraTrackTerm_anox_alpha~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_alpha_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_alpha_anox_relError),1,doseMatrixBixel.intraTrackTerm_alpha_anox_relError(doseMatrixBixel.intraTrackTerm_alpha_anox_relError~=0),dij.numOfVoxels,1);
   
    %% Electron
    dij.doseMatrixBixel(counterDijColumns).electronDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_electron),1,doseMatrixBixel.physicalDose_electron(doseMatrixBixel.physicalDose_electron~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).electronDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_electron_relError),1,doseMatrixBixel.physicalDose_electron_relError(doseMatrixBixel.physicalDose_electron_relError~=0),dij.numOfVoxels,1);
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).electronDosexRBE_aero{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_electron_aero),1,doseMatrixBixel.RBExDose_electron_aero(doseMatrixBixel.RBExDose_electron_aero~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).RBExDose_electron_aero_relError{1} = ...
        sparse(find(doseMatrixBixel.RBExDose_electron_aero_relError),1,doseMatrixBixel.RBExDose_electron_aero_relError(doseMatrixBixel.RBExDose_electron_aero_relError~=0),dij.numOfVoxels,1);
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_electron{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_electron),1,doseMatrixBixel.intraTrackTerm_aero_electron(doseMatrixBixel.intraTrackTerm_aero_electron~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_electron_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_electron_aero_relError),1,doseMatrixBixel.intraTrackTerm_electron_aero_relError(doseMatrixBixel.intraTrackTerm_electron_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).electronDosexRBE_anox{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_electron_anox),1,doseMatrixBixel.RBExDose_electron_anox(doseMatrixBixel.RBExDose_electron_anox~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_electron_anox_relError{1} = ...
%         sparse(find(doseMatrixBixel.RBExDose_electron_anox_relError),1,doseMatrixBixel.RBExDose_electron_anox_relError(doseMatrixBixel.RBExDose_electron_anox_relError~=0),dij.numOfVoxels,1);
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_electron{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_electron),1,doseMatrixBixel.intraTrackTerm_anox_electron(doseMatrixBixel.intraTrackTerm_anox_electron~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_electron_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_electron_anox_relError),1,doseMatrixBixel.intraTrackTerm_electron_anox_relError(doseMatrixBixel.intraTrackTerm_electron_anox_relError~=0),dij.numOfVoxels,1);
   
    %% Heavy ion
    dij.doseMatrixBixel(counterDijColumns).heavyIonDose{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_heavyIon),1,doseMatrixBixel.physicalDose_heavyIon(doseMatrixBixel.physicalDose_heavyIon~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).heavyIonDose_relError{1} = ...
        sparse(find(doseMatrixBixel.physicalDose_heavyIon_relError),1,doseMatrixBixel.physicalDose_heavyIon_relError(doseMatrixBixel.physicalDose_heavyIon_relError~=0),dij.numOfVoxels,1);
    % Check heavy ion RBE values
    if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
        if isfield(pln.propOpt, 'heavyIonRBEval_aero')
            RBEvalue_heavyIon_aero = pln.propOpt.heavyIonRBEval_aero;
        else
            RBEvalue_heavyIon_aero = 3.41;
        end
%         if isfield(pln.propOpt, 'heavyIonRBEval_anox')
%             RBEvalue_heavyIon_anox = pln.propOpt.heavyIonRBEval_anox;
%         else
%             RBEvalue_heavyIon_anox = 9.93;
%         end
    end
    % Aerobic environment
    dij.doseMatrixBixel(counterDijColumns).heavyIonDosexRBE_aero{1} = ...
        dij.doseMatrixBixel(counterDijColumns).heavyIonDose{1}*RBEvalue_heavyIon_aero;
    dij.doseMatrixBixel(counterDijColumns).RBExDose_heavyIon_aero_relError{1} = ...
        dij.doseMatrixBixel(counterDijColumns).heavyIonDose_relError{1};
    
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_heavyIon{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_aero_heavyIon),1,doseMatrixBixel.intraTrackTerm_aero_heavyIon(doseMatrixBixel.intraTrackTerm_aero_heavyIon~=0),dij.numOfVoxels,1);
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_heavyIon_relError{1} = ...
        sparse(find(doseMatrixBixel.intraTrackTerm_heavyIon_aero_relError),1,doseMatrixBixel.intraTrackTerm_heavyIon_aero_relError(doseMatrixBixel.intraTrackTerm_heavyIon_aero_relError~=0),dij.numOfVoxels,1);
%     % Anoxic environment
%     dij.doseMatrixBixel(counterDijColumns).heavyIonDosexRBE_anox{1} = ...
%         dij.doseMatrixBixel(counterDijColumns).heavyIonDose{1}*RBEvalue_heavyIon_anox;
%     dij.doseMatrixBixel(counterDijColumns).RBExDose_heavyIon_anox_relError{1} = ...
%         dij.doseMatrixBixel(counterDijColumns).heavyIonDose_relError{1};
%     
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_heavyIon{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_anox_heavyIon),1,doseMatrixBixel.intraTrackTerm_anox_heavyIon(doseMatrixBixel.intraTrackTerm_anox_heavyIon~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_heavyIon_relError{1} = ...
%         sparse(find(doseMatrixBixel.intraTrackTerm_heavyIon_anox_relError),1,doseMatrixBixel.intraTrackTerm_heavyIon_anox_relError(doseMatrixBixel.intraTrackTerm_heavyIon_anox_relError~=0),dij.numOfVoxels,1);
    
%     %% Photon dose
%     dij.doseMatrixBixel(counterDijColumns).photonDose{1} = ...
%         sparse(find(doseMatrixBixel.physicalDose_photonHeating),1,doseMatrixBixel.physicalDose_photonHeating(doseMatrixBixel.physicalDose_photonHeating~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).photonDose_relError{1} = ...
%         sparse(find(doseMatrixBixel.physicalDose_photonHeating_relError),1,doseMatrixBixel.physicalDose_photonHeating_relError(doseMatrixBixel.physicalDose_photonHeating_relError~=0),dij.numOfVoxels,1);
%     %% Neutron dose
%     dij.doseMatrixBixel(counterDijColumns).neutronDose{1} = ...
%         sparse(find(doseMatrixBixel.physicalDose_neutronHeating),1,doseMatrixBixel.physicalDose_neutronHeating(doseMatrixBixel.physicalDose_neutronHeating~=0),dij.numOfVoxels,1);
%     dij.doseMatrixBixel(counterDijColumns).neutronDose_relError{1} = ...
%         sparse(find(doseMatrixBixel.physicalDose_neutronHeating_relError),1,doseMatrixBixel.physicalDose_neutronHeating_relError(doseMatrixBixel.physicalDose_neutronHeating_relError~=0),dij.numOfVoxels,1);
end

%% B: Calculate alpha and beta values for neutron radiation using RMF model

% Calculate total dose
disp('*****')
disp('Total dose for calculation of RMF model parameters: sum of secondary charged particles only.')
disp('*****')
containerSummedDose = dij.doseMatrixBixel(counterDijColumns).protonDose{1} + ...
    dij.doseMatrixBixel(counterDijColumns).deuteronDose{1} + dij.doseMatrixBixel(counterDijColumns).tritonDose{1} + dij.doseMatrixBixel(counterDijColumns).alphaDose{1} + dij.doseMatrixBixel(counterDijColumns).he3Dose{1} + ...
    dij.doseMatrixBixel(counterDijColumns).heavyIonDose{1} + dij.doseMatrixBixel(counterDijColumns).electronDose{1};


%% B.1: Calculate dose weighted RBE for DSB induction
disp('*****')
disp('RBE weighted dose for calculation of RMF model parameters: sum of secondary charged particle dose x RBE(aero).')
disp('*****')
containerRBExDose_aero = dij.doseMatrixBixel(counterDijColumns).protonDosexRBE_aero{1} + ...
    dij.doseMatrixBixel(counterDijColumns).deuteronDosexRBE_aero{1} + dij.doseMatrixBixel(counterDijColumns).tritonDosexRBE_aero{1} + dij.doseMatrixBixel(counterDijColumns).alphaDosexRBE_aero{1} + dij.doseMatrixBixel(counterDijColumns).he3DosexRBE_aero{1} + ...
    dij.doseMatrixBixel(counterDijColumns).heavyIonDosexRBE_aero{1} + dij.doseMatrixBixel(counterDijColumns).electronDosexRBE_aero{1};

disp('*****')
disp('Dose weighted RBE values: ratio of RBExDose/Dose for each voxel (aerobic environment).')
disp('Secondary charged particles only.')
disp('*****')
containerDoseWeightedRBEvalues_aero = containerRBExDose_aero(containerSummedDose~=0)./containerSummedDose(containerSummedDose~=0);

% disp('*****')
% disp('RBE weighted dose for calculation of RMF model parameters: sum of secondary charged particle dose x RBE(anox).')
% disp('*****')
% containerRBExDose_anox = dij.doseMatrixBixel(counterDijColumns).protonDosexRBE_anox{1} + ...
%     dij.doseMatrixBixel(counterDijColumns).deuteronDosexRBE_anox{1} + dij.doseMatrixBixel(counterDijColumns).tritonDosexRBE_anox{1} + dij.doseMatrixBixel(counterDijColumns).alphaDosexRBE_anox{1} + dij.doseMatrixBixel(counterDijColumns).he3DosexRBE_anox{1} + ...
%     dij.doseMatrixBixel(counterDijColumns).heavyIonDosexRBE_anox{1} + dij.doseMatrixBixel(counterDijColumns).electronDosexRBE_anox{1};

% disp('*****')
% disp('Dose weighted RBE values: ratio of RBExDose/Dose for each voxel (anoxic environment).')
% disp('Secondary charged particles only.')
% disp('*****')
% containerDoseWeightedRBEvalues_anox = containerRBExDose_anox(containerSummedDose~=0)./containerSummedDose(containerSummedDose~=0);

%% B.2: Calculate dose weighted intra track term for RMF model
disp('*****')
disp('Dose weighted RBE values: ratio of zFxRBExRBExDose/Dose for each voxel.')
disp('Secondary charged particles only.')
disp('*****')

% Calculation of intra track term for all particles
containerIntraTrackTerm_aero = dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_proton{1} + ...
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_deuteron{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_triton{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_alpha{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_he3{1} + ...
    dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_heavyIon{1} ...
    + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_aero_electron{1};
% Calculation of intra track term for all particles
% containerIntraTrackTerm_anox = dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_proton{1} + ...
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_deuteron{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_triton{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_alpha{1} + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_he3{1} + ...
%     dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_heavyIon{1} ...
%     + dij.doseMatrixBixel(counterDijColumns).intraTrackTerm_anox_electron{1};


% Rescaling of zF intra track term according to zF rescaling factor
if isfield(pln, 'propMCDS')
    if isfield(pln.propMCDS, 'zFrescaleFactor')
        zFrescaleFactor = pln.propMCDS.zFrescaleFactor;
    end
else
        zFrescaleFactor = 1;
end

containerIntraTrackTerm_aero = containerIntraTrackTerm_aero * zFrescaleFactor;
% containerIntraTrackTerm_anox = containerIntraTrackTerm_anox * zFrescaleFactor;

containerDoseWeightedintraTrackTerm_aero =  spalloc(prod(ct.cubeDim), 1, prod(ct.cubeDim));
% containerDoseWeightedintraTrackTerm_anox =  spalloc(prod(ct.cubeDim), 1, prod(ct.cubeDim));
containerDoseWeightedintraTrackTerm_aero(containerSummedDose~=0) = containerIntraTrackTerm_aero(containerSummedDose~=0)./containerSummedDose(containerSummedDose~=0);
% containerDoseWeightedintraTrackTerm_anox(containerSummedDose~=0) = containerIntraTrackTerm_anox(containerSummedDose~=0)./containerSummedDose(containerSummedDose~=0);

% Check indices
if ~isempty(find(find(containerSummedDose) - find(containerDoseWeightedintraTrackTerm_aero)))
    error('Voxel indices of non-zero voxel elements of summed dose and intra track terms are not the same!')
% elseif ~isempty(find(find(containerSummedDose) - find(containerDoseWeightedintraTrackTerm_anox)))
%     error('Voxel indices of non-zero voxel elements of summed dose and intra track terms are not the same!')
elseif ~isempty(find(find(containerSummedDose) - find(containerRBExDose_aero)))
    error('Voxel indices of non-zero voxel elements of summed dose and RBE weighted dose are not the same!')
% elseif ~isempty(find(find(containerSummedDose) - find(containerRBExDose_anox)))
%     error('Voxel indices of non-zero voxel elements of summed dose and RBE weighted dose are not the same!')
end

%% B.3: Get alpha over beta ratio map for reference radiation
alphaBetaRatioMap = zeros(ct.cubeDim);

for counterRTstructures = 1:size(cst,1)
    % Check for predefined values
    if ~isfield(cst{counterRTstructures,5},'alphaBetaRatioX')
        cst{counterRTstructures,5}.alphaBetaRatioX = pln.propOpt.defaultLQmodel.ratioAlphaBeta;
    end
    if ~isfield(cst{counterRTstructures,5},'LQmodelPriority')
        if strcmpi(cst{counterRTstructures,2}, pln.propMCNP.bodyStructureName)
            cst{counterRTstructures,5}.LQmodelPriority = 0;
        elseif strcmp(cst{counterRTstructures,3}, 'OAR') && ~strcmpi(cst{counterRTstructures,2}, pln.propMCNP.bodyStructureName)
            cst{counterRTstructures,5}.LQmodelPriority = 1;
        elseif strcmp(cst{counterRTstructures,3}, 'TARGET')
            cst{counterRTstructures,5}.LQmodelPriority = 2;
        end
    end
end

for counterRTstructures = 1:size(cst,1)
    % Set values in alpha over beta ratio map for each voxel
    for counterOtherStruct = 1:size(cst,1)
        if size(unique([cst{counterRTstructures,4}{1}' cst{counterOtherStruct,4}{1}']),2)~= size(cst{counterRTstructures,4}{1},1)+size(cst{counterOtherStruct,4}{1},1) && counterOtherStruct~=counterRTstructures
            if cst{counterRTstructures,5}.LQmodelPriority > cst{counterOtherStruct,5}.LQmodelPriority
                alphaBetaRatioMap(cst{counterRTstructures,4}{1}) = cst{counterRTstructures,5}.alphaBetaRatioX;
            elseif cst{counterRTstructures,5}.LQmodelPriority < cst{counterOtherStruct,5}.LQmodelPriority
                alphaBetaRatioMap(cst{counterOtherStruct,4}{1}) = cst{counterOtherStruct,5}.alphaBetaRatioX;
            elseif cst{counterRTstructures,5}.LQmodelPriority == cst{counterOtherStruct,5}.LQmodelPriority
                if cst{counterRTstructures,5}.alphaBetaRatioX ~= cst{counterOtherStruct,5}.alphaBetaRatioX
                    warning('Tissue with identical LQ parameter priority but different alpha over beta values detected. Better re-define.')
                    alphaBetaRatioMap(cst{counterRTstructures,4}{1}) = cst{counterRTstructures,5}.alphaBetaRatioX;
                else
                    alphaBetaRatioMap(cst{counterRTstructures,4}{1}) = cst{counterRTstructures,5}.alphaBetaRatioX;
                end
            end
        else
            alphaBetaRatioMap(cst{counterRTstructures,4}{1}) = cst{counterRTstructures,5}.alphaBetaRatioX;
        end
    end
end

% Remember alpha/beta ratio map
dij.alphaBetaRatioMap{counterDijColumns} = alphaBetaRatioMap;

%% B.4: Calculate RBE_LD and RBE_HD values for neutron ray according to RMF model
dij.RBE_LD_aero{1}(:,counterDijColumns) = sparse(size(dij.physicalDose{1}(:,counterDijColumns),1),1);
dij.RBE_HD_aero{1}(:,counterDijColumns) = sparse(size(dij.physicalDose{1}(:,counterDijColumns),1),1);
% dij.RBE_LD_anox{1}(:,counterDijColumns) = sparse(size(dij.physicalDose{1}(:,counterDijColumns),1),1);
% dij.RBE_HD_anox{1}(:,counterDijColumns) = sparse(size(dij.physicalDose{1}(:,counterDijColumns),1),1);

% RBE_LD values
dij.RBE_LD_aero{1}(containerSummedDose~=0,counterDijColumns) = (containerDoseWeightedRBEvalues_aero+...
    (2./alphaBetaRatioMap(containerSummedDose~=0)).*containerDoseWeightedintraTrackTerm_aero(containerSummedDose~=0));

% RBE_HD values
dij.RBE_HD_aero{1}(containerSummedDose~=0,counterDijColumns) = containerDoseWeightedRBEvalues_aero;

% % RBE_LD values
% dij.RBE_LD_anox{1}(containerSummedDose~=0,counterDijColumns) = (containerDoseWeightedRBEvalues_anox+...
%     (2./alphaBetaRatioMap(containerSummedDose~=0)).*containerDoseWeightedintraTrackTerm_anox(containerSummedDose~=0));
% 
% % RBE_HD values
% dij.RBE_HD_anox{1}(containerSummedDose~=0,counterDijColumns) = containerDoseWeightedRBEvalues_anox;

end