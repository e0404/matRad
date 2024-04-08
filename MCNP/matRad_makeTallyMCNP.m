function matRad_makeTallyMCNP(ct, pln, fileID_C_rest, binIntervals)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description goes here.
%
% call
%   matRad_makeTallyMCNP(ct, pln)
%
% input
%   stf:
%
% output:
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User’s Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch pln.propMCNP.tallySpecifier
    case 'TotalDose_TMESH'
        switch pln.radiationMode
            case 'neutrons'
                disp('*****')
                disp('Tally type: TMESH3...')
                disp('*****')
                meshTally.typeCard = 'TMESH\n';
                meshTally.geometry = 'RMESH3 %s\n';
                meshTally.corA = 'CORA3 %.4f %dI %.4f\n';
                meshTally.corB = 'CORB3 %.4f %dI %.4f\n';
                meshTally.corC = 'CORC3 %.4f %dI %.4f\n';
                pln.propMCNP.tallyKeyword= 'TOTAL';
                
                %Write to text file
                fprintf(fileID_C_rest, 'C ***************************************************************\n');
                fprintf(fileID_C_rest, 'C C: Heating tally (one tally located in each voxel of the CT-data)\n');
                fprintf(fileID_C_rest, 'C ***************************************************************\n');
                fprintf(fileID_C_rest, meshTally.typeCard);
                fprintf(fileID_C_rest, meshTally.geometry, pln.propMCNP.tallyKeyword);
                fprintf(fileID_C_rest, meshTally.corA, .5*ct.resolution.x_resized, (ct.cubeDim(2)-1), ct.cubeDim(2)*ct.resolution.x_resized+.5*ct.resolution.x_resized);    % Caution: MATLAB indexing
                fprintf(fileID_C_rest, meshTally.corB, .5*ct.resolution.y_resized, (ct.cubeDim(1)-1), ct.cubeDim(1)*ct.resolution.y_resized+.5*ct.resolution.y_resized);
                fprintf(fileID_C_rest, meshTally.corC, .5*ct.resolution.z_resized, (ct.cubeDim(3)-1), ct.cubeDim(3)*ct.resolution.z_resized+.5*ct.resolution.z_resized);
                
                fprintf(fileID_C_rest, 'ENDMD\n');
                
                % Add tallies for RBE weighted dose
                if isfield(pln.propOpt,'bioOptimization')
                    if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
                        matRad_makeTallyMCNP4biolOptiRMF(ct, pln, fileID_C_rest, binIntervals)
                    else
                        disp('*****')
                        disp('Biological optimization mode unknown. Calculation will continue without biological weighting factors!')
                        disp('*****')
                    end
                end
            case 'neutrons_PLUS_photons'
                disp('*****')
                disp('Tally type: TMESH3...')
                disp('*****')
                meshTally.typeCard = 'TMESH\n';
                meshTally.geometry = 'RMESH3 %s\n';
                meshTally.corA = 'CORA3 %.4f %dI %.4f\n';
                meshTally.corB = 'CORB3 %.4f %dI %.4f\n';
                meshTally.corC = 'CORC3 %.4f %dI %.4f\n';
                pln.propMCNP.tallyKeyword= 'TOTAL';
                
                %Write to text file
                fprintf(fileID_C_rest, 'C ***************************************************************\n');
                fprintf(fileID_C_rest, 'C C: Heating tally (one tally located in each voxel of the CT-data)\n');
                fprintf(fileID_C_rest, 'C ***************************************************************\n');
                fprintf(fileID_C_rest, meshTally.typeCard);
                fprintf(fileID_C_rest, meshTally.geometry, pln.propMCNP.tallyKeyword);
                fprintf(fileID_C_rest, meshTally.corA, -.5*ct.resolution.x_resized, (ct.cubeDim(2)-1), ct.cubeDim(2)*ct.resolution.x_resized-.5*ct.resolution.x_resized);    % Caution: MATLAB indexing
                fprintf(fileID_C_rest, meshTally.corB, -.5*ct.resolution.y_resized, (ct.cubeDim(1)-1), ct.cubeDim(1)*ct.resolution.y_resized-.5*ct.resolution.y_resized);
                fprintf(fileID_C_rest, meshTally.corC, -.5*ct.resolution.z_resized, (ct.cubeDim(3)-1), ct.cubeDim(3)*ct.resolution.z_resized-.5*ct.resolution.z_resized);
                
                fprintf(fileID_C_rest, 'ENDMD\n');
                
                % Add tallies for RBE weighted dose
                if isfield(pln.propOpt,'bioOptimization')
                    if strcmp(pln.propOpt.bioOptimization,'RBExSecPartDose_MCDS_RMFmodel')
                        matRad_makeTallyMCNP4biolOptiRMF(ct, pln, fileID_C_rest, binIntervals)
                    else
                        disp('*****')
                        disp('Biological optimization mode unknown. Calculation will continue without biological weighting factors!')
                        disp('*****')
                    end
                end
        end
        
        fprintf(fileID_C_rest, 'PRINT 110\n');
        fprintf(fileID_C_rest, ['PRDMP ',num2str(ceil(pln.propMCNP.numberParticles)),' ',num2str(ceil(pln.propMCNP.numberParticles)), ' 1 ', num2str(ceil(pln.propMCNP.numberParticles)), '\n']);  % Control optional MCTAL textfile and set # dumps in RUNTPE to 1
        
        
        
        
    case 'KERMA_F4'
        switch pln.radiationMode
            case 'neutrons'
                % Get list of available tabulated KERMA factors
                kermaValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'KERMA_factors', filesep);
                kermaValues.kermaValuesList = dir([kermaValues.pathLocation, 'neutronKERMA*']); %
                
                if isfield(pln.propOpt,'bioOptimization')
                    switch pln.propOpt.bioOptimization
                        case 'const_RBExD_n'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_constFactor*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue(2,1) = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        case 'var_RBExD_n_ICRP'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_ICRP103*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        case 'var_RBExD_n_MCDS'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_StewardEtAl*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        otherwise
                            disp('*****')
                            disp('Biological optimization mode unknown. Calculation will continue without biological weighting factors!')
                            disp('*****')
                    end
                end
                
                for counterTally=1:size(binIntervals,2)-1
                    % Read KERMA Factors
                    fid_kermaVal = fopen([kermaValues.kermaValuesList(counterTally).folder, filesep, kermaValues.kermaValuesList(counterTally).name], 'r');
                    kermaValue = fscanf(fid_kermaVal, '%f', [2,inf]);
                    kermaValue(2,:) = kermaValue(2,:)*pln.propMCNP.normalizationFactor;   % multiply by factor s.th. no problems with double precision occur in MCNP
                    fclose(fid_kermaVal);
                    
                    %Prepare Writing Tally Cards
                    meshTally.typeCard = ['FMESH', int2str(counterTally), '014:N\n'];
                    meshTally.geometry = '        GEOM=xyz\n';
                    meshTally.origin = '        ORIGIN= %.13f %.13f %.13f\n';
                    meshTally.imesh = '        IMESH=%.13f\n';
                    meshTally.iints = '        IINTS=%d\n';
                    meshTally.jmesh = '        JMESH=%.13f\n';
                    meshTally.jints = '        JINTS=%d\n';
                    meshTally.kmesh = '        KMESH=%.13f\n';
                    meshTally.kints = '        KINTS=%d\n';
                    meshTally.deCard = ['DE', int2str(counterTally), '014 LIN\n'];
                    meshTally.deCard_input = '        %.13f\n';
                    meshTally.dfCard = ['DF', int2str(counterTally), '014 LIN\n'];
                    meshTally.dfCard_input = '        %.15f\n';
                    meshTally.commentCard = ['FC', int2str(counterTally), '014 F4 Tally recording Neutron KERMA by using Flux-to-Dose-Factors.\n'];
                    
                    disp('*****')
                    disp(['Tally type: FMESH4 with DE and DF card for ', pln.propMCNP.tallySpecifier, '. Output in [Gy] normalized to one particle.'])
                    disp('*****')
                    
                    % Write Tally Cards to Text File
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, 'C C.3: Tally (one tally located in each voxel of the CT-data)\n');
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, meshTally.typeCard);
                    fprintf(fileID_C_rest, meshTally.geometry);
                    fprintf(fileID_C_rest, meshTally.origin, -(ct.resolution.x_resized)/2, ...
                        -(ct.resolution.y_resized)/2,...
                        -(ct.resolution.z_resized)/2);
                    fprintf(fileID_C_rest, meshTally.imesh, (ct.cubeDim(2)*ct.resolution.x_resized)-ct.resolution.x_resized/2);
                    fprintf(fileID_C_rest, meshTally.iints, ct.cubeDim(2));
                    fprintf(fileID_C_rest, meshTally.jmesh, (ct.cubeDim(1)*ct.resolution.y_resized)-ct.resolution.y_resized/2);
                    fprintf(fileID_C_rest, meshTally.jints, ct.cubeDim(1));
                    fprintf(fileID_C_rest, meshTally.kmesh, (ct.cubeDim(3)*ct.resolution.z_resized)-ct.resolution.z_resized/2);
                    fprintf(fileID_C_rest, meshTally.kints, ct.cubeDim(3));
                    fprintf(fileID_C_rest, meshTally.deCard);
                    fprintf(fileID_C_rest, meshTally.deCard_input, kermaValue(1,:));
                    fprintf(fileID_C_rest, meshTally.dfCard);
                    if ~exist('RBEValue', 'var')  % Either use predefined RBE values or not
                        fprintf(fileID_C_rest, meshTally.dfCard_input, (kermaValue(2,:)*1e4)); % Multiplication by 1e4 since F4 tally detects flux in #/cm^2 and KERMA factors are given in Gy/m^2
                    elseif exist('RBEValue', 'var')
                        fprintf(fileID_C_rest, meshTally.dfCard_input, (RBEValue(2,:).*kermaValue(2,:)*1e4)); % Multiplication by 1e4 since F4 tally detects flux in #/cm^2 and KERMA factors are given in Gy/m^2
                    end
                    fprintf(fileID_C_rest, meshTally.commentCard);
                end
            case 'neutrons_PLUS_photons'
                % Get list of available tabulated KERMA factors
                kermaValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'KERMA_factors', filesep);
                kermaValues.kermaValuesList = dir([kermaValues.pathLocation, 'neutronKERMA*']); %
                
                if isfield(pln.propOpt,'bioOptimization')
                    switch pln.propOpt.bioOptimization
                        case 'const_RBExD_n'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_constFactor*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue(2,1) = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        case 'var_RBExD_n_ICRP'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_ICRP103*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        case 'var_RBExD_n_MCDS'
                            RBEValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'RBEfactors', filesep);
                            RBEValues.RBEValuesList = dir([RBEValues.pathLocation, 'neutronRBE_StewardEtAl*']);
                            fid_RBEVal = fopen([RBEValues.RBEValuesList.folder, filesep, RBEValues.RBEValuesList.name], 'r');
                            RBEValue = fscanf(fid_RBEVal, '%f', [2,inf]);
                            fclose(fid_RBEVal);
                        otherwise
                            disp('*****')
                            disp('Biological optimization mode unknown. Calculation will continue without biological weighting factors!')
                            disp('*****')
                    end
                end
                
                for counterTally=1:size(binIntervals,2)-1
                    % Read KERMA Factors
                    fid_kermaVal = fopen([kermaValues.kermaValuesList(counterTally).folder, filesep, kermaValues.kermaValuesList(counterTally).name], 'r');
                    kermaValue = fscanf(fid_kermaVal, '%f', [2,inf]);
                    kermaValue(2,:) = kermaValue(2,:)*pln.propMCNP.normalizationFactor;   % multiply by factor s.th. no problems with double precision occur in MCNP
                    fclose(fid_kermaVal);
                    
                    %Prepare Writing Tally Cards
                    meshTally.typeCard = ['FMESH', int2str(counterTally), '014:N\n'];
                    meshTally.geometry = '        GEOM=xyz\n';
                    meshTally.origin = '        ORIGIN= %.13f %.13f %.13f\n';
                    meshTally.imesh = '        IMESH=%.13f\n';
                    meshTally.iints = '        IINTS=%d\n';
                    meshTally.jmesh = '        JMESH=%.13f\n';
                    meshTally.jints = '        JINTS=%d\n';
                    meshTally.kmesh = '        KMESH=%.13f\n';
                    meshTally.kints = '        KINTS=%d\n';
                    meshTally.deCard = ['DE', int2str(counterTally), '014 LIN\n'];
                    meshTally.deCard_input = '        %.13f\n';
                    meshTally.dfCard = ['DF', int2str(counterTally), '014 LIN\n'];
                    meshTally.dfCard_input = '        %.15f\n';
                    meshTally.commentCard = ['FC', int2str(counterTally), '014 F4 Tally recording Neutron KERMA by using Flux-to-Dose-Factors.\n'];
                    
                    disp('*****')
                    disp(['Tally type: FMESH4 with DE and DF card for ', pln.propMCNP.tallySpecifier, '. Output in [Gy] normalized to one particle.'])
                    disp('*****')
                    
                    % Write Tally Cards to Text File
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, 'C C.3: Tally (one tally located in each voxel of the CT-data)\n');
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, meshTally.typeCard);
                    fprintf(fileID_C_rest, meshTally.geometry);
                    fprintf(fileID_C_rest, meshTally.origin, -(ct.resolution.x_resized)/2, ...
                        -(ct.resolution.y_resized)/2,...
                        -(ct.resolution.z_resized)/2);
                    fprintf(fileID_C_rest, meshTally.imesh, (ct.cubeDim(2)*ct.resolution.x_resized)-ct.resolution.x_resized/2);
                    fprintf(fileID_C_rest, meshTally.iints, ct.cubeDim(2));
                    fprintf(fileID_C_rest, meshTally.jmesh, (ct.cubeDim(1)*ct.resolution.y_resized)-ct.resolution.y_resized/2);
                    fprintf(fileID_C_rest, meshTally.jints, ct.cubeDim(1));
                    fprintf(fileID_C_rest, meshTally.kmesh, (ct.cubeDim(3)*ct.resolution.z_resized)-ct.resolution.z_resized/2);
                    fprintf(fileID_C_rest, meshTally.kints, ct.cubeDim(3));
                    fprintf(fileID_C_rest, meshTally.deCard);
                    fprintf(fileID_C_rest, meshTally.deCard_input, kermaValue(1,:));
                    fprintf(fileID_C_rest, meshTally.dfCard);
                    if ~exist('RBEValue', 'var')  % Either use predefined RBE values or not
                        fprintf(fileID_C_rest, meshTally.dfCard_input, (kermaValue(2,:)*1e4)); % Multiplication by 1e4 since F4 tally detects flux in #/cm^2 and KERMA factors are given in Gy/m^2
                    elseif exist('RBEValue', 'var')
                        fprintf(fileID_C_rest, meshTally.dfCard_input, (RBEValue(2,:).*kermaValue(2,:)*1e4)); % Multiplication by 1e4 since F4 tally detects flux in #/cm^2 and KERMA factors are given in Gy/m^2
                    end
                    fprintf(fileID_C_rest, meshTally.commentCard);
                end
                kermaValues.kermaValuesList = dir([kermaValues.pathLocation, 'photonKERMA*']); %
                for counterTally=1:size(binIntervals,2)-1
                    % Read KERMA Factors
                    fid_kermaVal = fopen([kermaValues.kermaValuesList(counterTally).folder, filesep, kermaValues.kermaValuesList(counterTally).name], 'r');
                    kermaValue = fscanf(fid_kermaVal, '%f', [2,inf]);
                    kermaValue(2,:) = kermaValue(2,:)*pln.propMCNP.normalizationFactor;   % multiply by factor s.th. no problems with double precision occur in MCNP
                    fclose(fid_kermaVal);
                    
                    %Prepare Writing Tally Cards
                    meshTally.typeCard = ['FMESH', int2str(counterTally), '024:P\n'];
                    meshTally.geometry = '        GEOM=xyz\n';
                    meshTally.origin = '        ORIGIN= %.13f %.13f %.13f\n';
                    meshTally.imesh = '        IMESH=%.13f\n';
                    meshTally.iints = '        IINTS=%d\n';
                    meshTally.jmesh = '        JMESH=%.13f\n';
                    meshTally.jints = '        JINTS=%d\n';
                    meshTally.kmesh = '        KMESH=%.13f\n';
                    meshTally.kints = '        KINTS=%d\n';
                    meshTally.deCard = ['DE', int2str(counterTally), '024 LIN\n'];
                    meshTally.deCard_input = '        %.13f\n';
                    meshTally.dfCard = ['DF', int2str(counterTally), '024 LIN\n'];
                    meshTally.dfCard_input = '        %.15f\n';
                    meshTally.commentCard = ['FC', int2str(counterTally), '024 F4 Tally recording Neutron KERMA by using Flux-to-Dose-Factors.\n'];
                    
                    disp('*****')
                    disp(['Tally type: FMESH4 with DE and DF card for ', pln.propMCNP.tallySpecifier, '. Output in [Gy] normalized to one particle.'])
                    disp('*****')
                    
                    % Write Tally Cards to Text File
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, 'C C.3: Tally (one tally located in each voxel of the CT-data)\n');
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, meshTally.typeCard);
                    fprintf(fileID_C_rest, meshTally.geometry);
                    fprintf(fileID_C_rest, meshTally.origin, -(ct.resolution.x_resized)/2, ...
                        -(ct.resolution.y_resized)/2,...
                        -(ct.resolution.z_resized)/2);
                    fprintf(fileID_C_rest, meshTally.imesh, (ct.cubeDim(2)*ct.resolution.x_resized)-ct.resolution.x_resized/2);
                    fprintf(fileID_C_rest, meshTally.iints, ct.cubeDim(2));
                    fprintf(fileID_C_rest, meshTally.jmesh, (ct.cubeDim(1)*ct.resolution.y_resized)-ct.resolution.y_resized/2);
                    fprintf(fileID_C_rest, meshTally.jints, ct.cubeDim(1));
                    fprintf(fileID_C_rest, meshTally.kmesh, (ct.cubeDim(3)*ct.resolution.z_resized)-ct.resolution.z_resized/2);
                    fprintf(fileID_C_rest, meshTally.kints, ct.cubeDim(3));
                    fprintf(fileID_C_rest, meshTally.deCard);
                    fprintf(fileID_C_rest, meshTally.deCard_input, kermaValue(1,:));
                    fprintf(fileID_C_rest, meshTally.dfCard);
                    fprintf(fileID_C_rest, meshTally.dfCard_input, kermaValue(2,:)*1e4);
                    fprintf(fileID_C_rest, meshTally.commentCard);
                end
            case 'photons'
                % Get list of available tabulated KERMA factors
                kermaValues.pathLocation = fullfile(matRad_getMATRADdirectory, 'MCNP', 'KERMA_factors', filesep);
                kermaValues.kermaValuesList = dir([kermaValues.pathLocation, 'photonKERMA*']); %
                for counterTally=1:size(binIntervals,2)-1
                    % Read KERMA Factors
                    fid_kermaVal = fopen([kermaValues.kermaValuesList(counterTally).folder, filesep, kermaValues.kermaValuesList(counterTally).name], 'r');
                    kermaValue = fscanf(fid_kermaVal, '%f', [2,inf]);
                    kermaValue(2,:) = kermaValue(2,:)*pln.propMCNP.normalizationFactor;   % multiply by factor s.th. no problems with double precision occur in MCNP
                    fclose(fid_kermaVal);
                    
                    %Prepare Writing Tally Cards
                    meshTally.typeCard = ['FMESH', int2str(counterTally), '024:P\n'];
                    meshTally.geometry = '        GEOM=xyz\n';
                    meshTally.origin = '        ORIGIN= %.13f %.13f %.13f\n';
                    meshTally.imesh = '        IMESH=%.13f\n';
                    meshTally.iints = '        IINTS=%d\n';
                    meshTally.jmesh = '        JMESH=%.13f\n';
                    meshTally.jints = '        JINTS=%d\n';
                    meshTally.kmesh = '        KMESH=%.13f\n';
                    meshTally.kints = '        KINTS=%d\n';
                    meshTally.deCard = ['DE', int2str(counterTally), '024 LIN\n'];
                    meshTally.deCard_input = '        %.13f\n';
                    meshTally.dfCard = ['DF', int2str(counterTally), '024 LIN\n'];
                    meshTally.dfCard_input = '        %.15f\n';
                    meshTally.commentCard = ['FC', int2str(counterTally), '024 F4 Tally recording Neutron KERMA by using Flux-to-Dose-Factors.\n'];
                    
                    disp('*****')
                    disp(['Tally type: FMESH4 with DE and DF card for ', pln.propMCNP.tallySpecifier, '. Output in [Gy] normalized to one particle.'])
                    disp('*****')
                    
                    % Write Tally Cards to Text File
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, 'C C.3: Tally (one tally located in each voxel of the CT-data)\n');
                    fprintf(fileID_C_rest, 'C ***************************************************************\n');
                    fprintf(fileID_C_rest, meshTally.typeCard);
                    fprintf(fileID_C_rest, meshTally.geometry);
                    fprintf(fileID_C_rest, meshTally.origin, -(ct.resolution.x_resized)/2, ...
                        -(ct.resolution.y_resized)/2,...
                        -(ct.resolution.z_resized)/2);
                    fprintf(fileID_C_rest, meshTally.imesh, (ct.cubeDim(2)*ct.resolution.x_resized)-ct.resolution.x_resized/2);
                    fprintf(fileID_C_rest, meshTally.iints, ct.cubeDim(2));
                    fprintf(fileID_C_rest, meshTally.jmesh, (ct.cubeDim(1)*ct.resolution.y_resized)-ct.resolution.y_resized/2);
                    fprintf(fileID_C_rest, meshTally.jints, ct.cubeDim(1));
                    fprintf(fileID_C_rest, meshTally.kmesh, (ct.cubeDim(3)*ct.resolution.z_resized)-ct.resolution.z_resized/2);
                    fprintf(fileID_C_rest, meshTally.kints, ct.cubeDim(3));
                    fprintf(fileID_C_rest, meshTally.deCard);
                    fprintf(fileID_C_rest, meshTally.deCard_input, kermaValue(1,:));
                    fprintf(fileID_C_rest, meshTally.dfCard);
                    fprintf(fileID_C_rest, meshTally.dfCard_input, kermaValue(2,:)*1e4);
                    fprintf(fileID_C_rest, meshTally.commentCard);
                end
        end
        fprintf(fileID_C_rest, 'PRINT 110\n');      % print table 110 to check starting particles
        fprintf(fileID_C_rest, 'PRDMP j j j 1\n');  % Set # dumps in RUNTPE to 1
        % fprintf(fileID_C_rest, 'SPDTL FORCE\n');    % force MCNP to run lattice speed tally enhancement
end