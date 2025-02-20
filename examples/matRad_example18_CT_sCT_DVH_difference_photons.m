%% Example: Photon Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% In this example, we will show:
% (i) how to load DICOM patient data into matRad and modify dosimetric
% objectives and constraints;
% (ii) how to set up a photon dose calculation and optimization based on the real CT;
% (iii) how to translate the optimized plan to be applied to a synthetic CT image; and
% (iv) how to calculate the DVH differences between the plans.

%% Patient Data Import from DICOM  
% Start with a clean MATLAB environment, then import the real CT data of a patient.  
% For this purpose, LIVER phantom data has been adapted and saved as a DICOM file with a structure set.  
% Functions for importing DICOM files will search for .dcm files in the directory  
% and automatically recognize 3D volumes and structure files.  
% The DICOM files will then be read into the workspace as 'ct' and 'cst',  
% representing the CT images and the structure set, respectively.  
% Ensure that the matRad root directory, along with all its subdirectories,  
% is added to the MATLAB search path.  

clear;
matRad_cfg = matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
patDir= "matRad/phantoms/synthetic_ct/"; % If you want to export your data further, use "userdata" folder
name="pat001";
patDir_real_CT= fullfile(patDir,name,"real_CT");
imp_real_CT = matRad_DicomImporter(patDir_real_CT);
matRad_importDicom(imp_real_CT);
% Review the exported file and automatically identified  
% optimization objectives and constraints. 
matRadGUI;

%% Modifying Plan Optimization Objectives  
% The constraints for all defined volumes of interest (VOIs) are stored in the cst cell.  
% When importing DICOM structure files, matRad uses matRad_createCst() to generate this cell.  
% Default regular expressions are used to define target objectives. However,  
% additional constraints for organs at risk (OARs) need to be added for the plan.  
% This can be done by directly modifying the cst cell.  
% If different volumes with varying optimization objectives need to be loaded,  
% refer to the documentation for matRad_createCst() and choose suitable options  
% to create a custom implementation. In this case, the plan will be modified directly here.  

for i = 1:size(cst, 1)
    % Accessing each optimization objective and constain
    % Read more about how to set it in matRad_DoseOptimizationFunction - Superclass for objectives and constraints. 
    % This is the superclass for all objectives and constraints to enable easy 
    % identification.

    if cst{i,3} == "OAR"
        if ~isempty(regexpi(cst{i,2},'spinal', 'once'))
            objective = DoseObjectives.matRad_MaxDVH; 
            objective.penalty = 1;
            objective.parameters = {12,1};  %dose, to volume
            cst{i,6} = {};
            cst{i,6}{1} = struct(objective);

        elseif ~isempty(regexpi(cst{i,2},'stomach','once')) || ...
                ~isempty(regexpi(cst{i,2},'duodenum', 'once')) 
            objective = DoseObjectives.matRad_MaxDVH;
            objective.penalty = 1;
            objective.parameters = {30,1};  %dose, to volume
            cst{i,6} = {};
            cst{i,6}{1} = struct(objective);

        end
    end
end

% Verify that the new objectives have been added and are visible in the
% user interface. We save it for further synthetic CT calculations
matRadGUI;
real_CT_cst=cst;

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this case
% we want to use photons. Then, we need to define a treatment machine to 
% correctly load the corresponding base data. matRad includes base data for
% generic photon linear accelerator called 'Generic'. By this means matRad 
% will look for 'photons_Generic.mat' in our root directory and will use 
% the data provided in there for dose calculation.
% The number of fractions is set to 30. Internally, matRad considers the 
% fraction dose for optimization, however, objetives and constraints are 
% defined for the entire treatment.

pln.radiationMode   = 'photons';  
pln.machine         = 'Generic';
pln.numOfFractions  = 1;

%%
% Define the biological model used for modeling biological dose (esp. for
% particles).
% Possible biological models are:
% none:        use no specific biological model
% constRBE:    use a constant RBE
% MCN:         use the variable RBE McNamara model for protons
% WED:         use the variable RBE Wedenberg model for protons
% LEM:         use the biophysical variable RBE Local Effect model for carbons
% As we are  using photons, we simply set the parameter to 'none'
pln.bioModel = 'none';

%% 
% It is possible to request multiple error scenarios for robustness
% analysis and optimization. Here, we just use the "nominal scenario"
% (nomScen)
pln.multScen = 'nomScen';

%%
% Now we have to set some beam parameters. We can define multiple beam 
% angles for the treatment which might be used in generic fashion for quick dose calculation. All corresponding couch angles are set to 0 at this 
% point. Moreover, we set the bixelWidth to 4, which results in a beamlet 
% size of 4 x 4 mm in the isocenter plane. 

pln.propStf.gantryAngles   = [0 50 100 150 200 250 300];
pln.propStf.couchAngles    = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.bixelWidth     = 4;

% Obtain the number of beams and voxels from the existing variables and 
% calculate the iso-center which is per default the center of gravity of 
% all target voxels.
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = matRad_getIsoCenter(cst,ct,0);

%% Dose calculation settings
% set resolution of dose calculation and optimization
pln.propDoseCalc.doseGrid.resolution.x = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 4; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 4; % [mm]

%%
% Enable sequencing and disable direct aperture optimization (DAO) for now.
% A DAO optimization is shown in a seperate example.
pln.propSeq.runSequencing = 1;
pln.propOpt.runDAO        = 0;


%%
% and et voila our treatment plan structure is ready. Lets have a look:
%disp(pln);


%% Generate Beam Geometry STF
% The steering file struct comprises the complete beam geometry along with 
% ray position, pencil beam positions and energies, source to axis distance (SAD) etc.

stf = matRad_generateStf(ct,cst,pln);


%% Dose Calculation for real CT
% Let's generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities for real CT of a patient. Having dose influences available 
% allows subsequent inverse optimization.
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for IMRT for real CT
% The goal of the fluence optimization is to find a set of beamlet/pencil 
% beam weights which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger once the GUI to 
% visualize the optimized dose cubes.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

% Get the weights of the optimized plan to apply them to the synthetic CT image of the patient later.  
% Call the matRad_planAnalysis function with the prepared arguments to extract DVH  
% parameters calculated for the real CT image. 
weights=resultGUI.w;
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);

%Get the plan parameters
dvh = resultGUI.dvh;
qi = resultGUI.qi;

% Save DVH parameters calculated on the real CT to the table  
% for further comparison with plans calculated on the synthetic CT.  
% Include the patient number and indicate the CT type as "real"  
% to facilitate delta calculations between the plans later.  
dvh_table_real=struct2table(qi);
% Select only DVH parameters from QI table you are interested in comparison
dvh_table_real=dvh_table_real(:,1:9);
dvh_table_real.patient= repmat(char(name),length(qi),1);
dvh_table_real.ct_type = repmat(char("real"),length(qi),1);
% Check DVH table for real CT
disp(dvh_table_real);

%% Now clear the data from the real CT image, except for the plan parameters,  
% and load the synthetic (fake) CT image of the same patient.  
% It is important that your image sets are compatible (i.e., same number of CT slices,  
% same isocenter position, etc.). We will re-use the structure file from real CT with adjusted objectives 
% for dose calculations.
clear resultGUI ct cst idx qi* dvh dij;
patDir_fake_CT= fullfile(patDir,name,"fake_CT");
imp_fake_CT = matRad_DicomImporter(patDir_fake_CT);
matRad_importDicom(imp_fake_CT);
delete(matRadGUI);
cst=real_CT_cst;
% Review the exported file and previously identified
% optimization objectives and constraints. 
matRadGUI;

%% Perform the dose calculation for synthetic CT by using the weights of plan, calculated on real CT 
% (i.e., obtain the dij variable) for the synthetic CT image. 
%TBD ask a question whether it is correct in the new version
resultGUI = matRad_calcDoseDirect(ct,stf,pln,cst, weights);
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);
matRadGUI;
%Get the plan parameters calculated on synthetic CT
dvh = resultGUI.dvh;
qi = resultGUI.qi;

% In the similar fashion, save DVH parameters calculated on the synthetic CT to the table  
dvh_table_fake=struct2table(qi);
% Select the same DVH parameters for further comparison
dvh_table_fake=dvh_table_fake(:,1:9);
dvh_table_fake.patient= repmat(char(name),length(qi),1);
dvh_table_fake.ct_type = repmat(char("fake"),length(qi),1);

%% Calculate the difference in between of DVH calculated on real CT (dvh_table_real) and synthetic CT (dvh_table_fake)
dvh_table_diff=dvh_table_real;
for i =1:height(dvh_table_real)
    for j=["mean","std","max", "min", "D_2","D_5", "D_95", "D_98"]
        if cell2mat(dvh_table_real{i,"name"})==cell2mat(dvh_table_fake{i,"name"})
            dvh_table_diff{i,j}=(dvh_table_real{i,j}-dvh_table_fake{i,j})/dvh_table_real{i,j}*100;
        else
            dvh_table_diff{i,j}="error";
        end
    end
end

%% Check the difference
disp(dvh_table_diff)

%%% Save the results to CSV file
%file_path_dvh_diff = "YOUR PATH"
%writetable(dvh_table_diff, file_path_dvh_diff);

close all;
clear;