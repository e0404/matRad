function [control_makeTargetMCNP, fileID_A, fileID_B, geometryOption] = matRad_makeTargetMCNP(ct, simPropMCNP)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write target as MCNP input file using ct data, predefined simulation
% properties and tissue specific HU intervals defined for segmentation
%
% call
%   [control_makeTargetMCNP, fileID_A, fileID_B] = matRad_makeTargetMCNP(ct, simPropMCNP)
%
% input
%   ct:             ct cube
%   simPropMCNP:    includes .MCNP_limitNumberOfElements and .loopCounter
%                   to control number of elements and loop counter to
%                   rerun creation if not successful
%
% output
%   control_makeTargetMCNP: Control variable to check if target was created
%   sucessfully
%   fileID_A/B:     File IDs to control block A and B txt-files
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User’s Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg =  MatRad_Config.instance();

%% Check if geometry option is defined, otherwise set to default, i.e. use repeated structure option (lattice)
if ~isfield(simPropMCNP, 'geometryOption')
    simPropMNCP.geometryOption = 'Lattice';
    geometryOption =  'Lattice'; % dummy to give back by function
end

%% Define voxel phantom for MCNP calculation
switch simPropMNCP.geometryOption
    case 'Lattice'
        %% Generate universe matrix and open text files
        universeMatrix = zeros(ct.doseGridCT.cubeDim, 'uint8');    % create a 3D matrix and fill with material/universe indices
        % note: every universe corresponds to one material
        for counterMaterial = 1:size(ct.doseGridCT.tissueBin,2)
            universeMatrix(ct.doseGridCT.tissueBin(counterMaterial).linIndVol) = counterMaterial;  % fill
        end
        clear counterMaterial
        
        % Generate log-file    
        pathRunfiles = strcat(matRad_cfg.matRadRoot,filesep,'MCNP', filesep, 'runfiles_tmp', filesep);
        fileID_A = fopen(strcat(pathRunfiles,'blockA.txt'), 'w');
        fileID_B = fopen(strcat(pathRunfiles,'blockB.txt'), 'w');
        
        %% Write block A: Cells
        fprintf(fileID_A, 'MCNP Runfile for Dose Calculation\n');
        fprintf(fileID_A, 'C ***************************************************************\n');
        fprintf(fileID_A, 'C ***************************************************************\n');
        fprintf(fileID_A, 'C Block A: Cells\n');
        fprintf(fileID_A, 'C ***************************************************************\n');
        
        % Start with definition of cells where cell index equals material index equals universe index
        blockA.universeSpec = '%d %d %.9f %d u=%d $ Material %d\n';
        blockA.universeSpecLikeBut = '%d like 1 but mat=%d rho=%.9f u=%d $ Material %d\n';
        for cellCounter = 1:size(ct.doseGridCT.tissueBin,2)
            if ~isnan(mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol)))
                if cellCounter==1
                    fprintf(fileID_A, blockA.universeSpec, ...
                        cellCounter, ...        % cell number
                        cellCounter, ...        % material number
                        -mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol)), ...    % material density
                        -(100+cellCounter), ... % surface
                        cellCounter, ...        % universe number
                        cellCounter);           % material number for comment in runfile
                else
                    fprintf(fileID_A, blockA.universeSpecLikeBut, ...
                        cellCounter, ...        % cell number
                        cellCounter, ...        % material number
                        -mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol)), ...    % material density
                        cellCounter, ...        % universe number
                        cellCounter);           % material number for comment in runfile
                end
            elseif isnan(mean(ct.doseGridCT.density{1,1}(ct.doseGridCT.tissueBin(cellCounter).linIndVol)))
                if cellCounter==1
                    fprintf(fileID_A, blockA.universeSpec, ...
                        cellCounter, ...        % cell number
                        cellCounter, ...        % material number
                        1,...                   % set material density to one
                        -(100+cellCounter), ... % surface
                        cellCounter, ...        % universe number
                        cellCounter);           % material number for comment in runfile
                else
                    fprintf(fileID_A, blockA.universeSpecLikeBut, ...
                        cellCounter, ...        % cell number
                        cellCounter, ...        % material number
                        1,...                   % set material density to one
                        cellCounter, ...        % universe number
                        cellCounter);           % material number for comment in runfile
                end
            end
        end
        
        % Define lattice cell and fill with universes
        cellCounter = cellCounter +1;
        blockA.latticeSpec = '%d 0 -505 506 -507 508 -509 510 lat=1 u=%d fill=0:%d 0:%d 0:%d\n'; % surfaces 505:510 always define lattice cell [0,0,0]
        fprintf(fileID_A, blockA.latticeSpec, ...
            cellCounter, ...
            cellCounter, ...
            ct.doseGridCT.cubeDim(2)-1, ...
            ct.doseGridCT.cubeDim(1)-1, ...
            ct.doseGridCT.cubeDim(3)-1);
        for i3 = 1:ct.doseGridCT.cubeDim(3) % write universe index for each cell
            for i1 = 1:ct.doseGridCT.cubeDim(1)
                for i2 = 1:ct.doseGridCT.cubeDim(2)
                    fprintf(fileID_A, ['           ', int2str(universeMatrix(i1,i2,i3)), '\n']);
                end
            end
        end

        % Define inner world for universe to position ct grid
        cellCounter = cellCounter +1;
        blockA.innerWorldSpec = '%d 0 -911 fill=%d\n';
        fprintf(fileID_A, blockA.innerWorldSpec, cellCounter, cellCounter-1);

        % Define inner world layer around ct grid to simulate additional air
        cellCounter = cellCounter +1;
        blockA.innerWorldSpec = '%d 1 -0.001205000 911 -912\n';
        fprintf(fileID_A, blockA.innerWorldSpec, cellCounter);


        % Define outer world as particle graveyard
        cellCounter = cellCounter +1;
        blockA.outerWorldSpec = '%d 0 912\n';
        fprintf(fileID_A, blockA.outerWorldSpec, cellCounter);
        
        clear cellCounter
        
        fprintf(fileID_A, '\n');
        fclose(fileID_A);
        
        %% Write block B: Surfaces
        fprintf(fileID_B, 'C ***************************************************************\n');
        fprintf(fileID_B, 'C ***************************************************************\n');
        fprintf(fileID_B, 'C Block B: Surfaces\n');
        fprintf(fileID_B, 'C ***************************************************************\n');
        
        % Write universe voxel: cell 101
        blockB.universeVoxelSpec = '101 RPP %.13f %.13f %.13f %.13f %.13f %.13f $ voxel surface for universe 1\n';
        fprintf(fileID_B, blockB.universeVoxelSpec, ...
                0.25*ct.doseGridCT.x_MCNP, ...
                1.75*ct.doseGridCT.x_MCNP, ...
                0.25*ct.doseGridCT.y_MCNP, ...
                1.75*ct.doseGridCT.y_MCNP, ...
                0.25*ct.doseGridCT.z_MCNP, ...
                1.75*ct.doseGridCT.z_MCNP);
        
        % Write surfaces for lattice orientation: lattice cell (0,0,0)
        fprintf(fileID_B, ['505 PX ', num2str(ct.doseGridCT.x_MCNP*1.5), ' $ defines x-direction\n']);
        fprintf(fileID_B, ['506 PX ', num2str(ct.doseGridCT.x_MCNP/2), '\n']);
        fprintf(fileID_B, ['507 PY ', num2str(ct.doseGridCT.y_MCNP*1.5), '$ defines y-direction\n']);
        fprintf(fileID_B, ['508 PY ', num2str(ct.doseGridCT.y_MCNP/2), '\n']);
        fprintf(fileID_B, ['509 PZ ', num2str(ct.doseGridCT.z_MCNP*1.5), '$ defines z-direction\n']);
        fprintf(fileID_B, ['510 PZ ', num2str(ct.doseGridCT.z_MCNP/2), '\n']);
        
        
        
        % Write additional surfaces for universe box: 
        latticeMargin = 0.005; % [cm] -> avoid problems with outer surfaces of lattice, otherwise particles might get lost
        fprintf(fileID_B, ['911 RPP ', ... 
            num2str(latticeMargin+ct.doseGridCT.x_MCNP/2), ' ',...
            num2str(ct.doseGridCT.cubeDim(2)*ct.doseGridCT.x_MCNP+ct.doseGridCT.x_MCNP/2-latticeMargin),' ',...
            num2str(latticeMargin+ct.doseGridCT.y_MCNP/2), ' ',...
            num2str(ct.doseGridCT.cubeDim(1)*ct.doseGridCT.y_MCNP+ct.doseGridCT.y_MCNP/2-latticeMargin), ' ',...
            num2str(latticeMargin+ct.doseGridCT.z_MCNP/2), ' ',...
            num2str(ct.doseGridCT.cubeDim(3)*ct.doseGridCT.z_MCNP+ct.doseGridCT.z_MCNP/2-latticeMargin),'\n']);

        % Write additional surfaces for layer around ct grid
        additionalLayerThickness = 180;         % this value defines the maximum source-to-axis (SAD) distance that can be simulated before particles start in imp:0 cell, 
                                                % holds only for 0degree couch rotation and gantry 0, 90, 180, 270 degree, shorter SAD necessary otherwise due to rotated surface source        
        
        fprintf(fileID_B, ['912 RPP ',num2str(-(ct.doseGridCT.x_MCNP/2+additionalLayerThickness)), ' ',...
             num2str(ct.doseGridCT.cubeDim(2)*ct.doseGridCT.x_MCNP-ct.doseGridCT.x_MCNP/2+additionalLayerThickness), ' ',...
             num2str(-(ct.doseGridCT.y_MCNP/2+additionalLayerThickness)), ' ',...
             num2str(ct.doseGridCT.cubeDim(1)*ct.doseGridCT.y_MCNP-ct.doseGridCT.y_MCNP/2+additionalLayerThickness), ' ', ...
             num2str(-(ct.doseGridCT.z_MCNP/2+additionalLayerThickness)), ' ', ...
             num2str(ct.doseGridCT.cubeDim(3)*ct.doseGridCT.z_MCNP-ct.doseGridCT.z_MCNP/2+additionalLayerThickness), '\n']);
                
        % Write dummy surface to position RSSA input:
        fprintf(fileID_B, '1001 1 PX 0 $ dummy surface for RSSA positioning\n'); 

                      
        fprintf(fileID_B, '\n');
        
        fclose(fileID_B);
        
        control_makeTargetMCNP =true;  
end