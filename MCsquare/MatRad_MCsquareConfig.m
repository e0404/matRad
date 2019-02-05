classdef MatRad_MCsquareConfig
% MatRad_MCsquareConfig class definition
% 
%
% References
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


    
    properties       
        
        %%% Simulation parameters:
        Num_Threads   =	0;		% Number of parallel calculation threads. Default: 0 = max available threads
        RNG_Seed      =	0;		% Seed for the random number generator (deterministic result only with single thread). Default: 0 = seed based on the time
        Num_Primaries = 1e6;		% Number of primary protons to simulate. Default: 1e7
        E_Cut_Pro     =	0.5;		% Energy cut (in MeV) below which heavy charged particles are locally absorbed. Default: 0.5
        D_Max	      =	0.2;		% Maximum distance between two step (cm). Default: 0.2
        Epsilon_Max   =	0.25;		% Fractional energy loss (dE/T) per step. Default: 0.25
        Te_Min	      =	0.05;		% Threshold energy (MeV) for the production of secondary electrons (currently locally absorbed). Default: 0.05
        % As a reference: 200 MeV protons can transfer a maximum energy of 0.5 MeV to ?-electrons which correspond to a range of 7 mm in lung tissues.
        
        %%% Input files
        CT_File                     = 'Patient.mhd';				% Name of the CT file. Default: CT.mhd
        HU_Density_Conversion_File	= 'Scanners/matRad_water/HU_Density_Conversion.txt';	% Name of the file containing HU to density conversion data. Default: HU_Density_Conversion.txt
        HU_Material_Conversion_File	= 'Scanners/matRad_water/HU_Material_Conversion.txt';	% Name of the file containing HU to material conversion data. Default: HU_Material_Conversion.txt
        BDL_Machine_Parameter_File  = 'BDL/BDL_matrad.txt';			% Name of the machine parameter file for the beam data library. Default: BDL.txt
        BDL_Plan_File               = 'PlanPencil.txt';			% Name of the plan file for the beam data library. Default: Plan.txt
        
        %%% Physical parameters
        Simulate_Nuclear_Interactions = true;     % Enable/Disable the simulation of nuclear interactions. Default: True
        Simulate_Secondary_Protons	  = true;         % Enable/Disable the simulation of secondary protons (emitted during nuclear interactions). Default: True
        Simulate_Secondary_Deuterons  = true;        % Enable/Disable the simulation of secondary deuterons (emitted during nuclear interactions). Default: True
        Simulate_Secondary_Alphas     = true;           % Enable/Disable the simulation of secondary alphas (emitted during nuclear interactions). Default: True
        
        
        %%% 4D simulation
        fourD_Mode			    = false;	% Enable/Disable the 4D simulation mode. Default: False
        fourD_Dose_Accumulation = false;		% Enable/Disable the dose accumulation for all 4D-CT phases. Default: False
        Field_type              = 'Velocity';	% Field type: Displacement or Velocity. Default: Velocity
        Create_Ref_from_4DCT	= false;		% Create the reference phase image from 4D CT images (True), or import the reference image (False). Default: False
        Create_4DCT_from_Ref	= false;		% Create 4D CT images by deforming the reference phase image (True), or import 4D CT images (False). Default: False
        Dynamic_delivery        = false;		% Enable/Disable simulation of dynamic delivery (interplay simulation). Default: False
        Breathing_period        = 7.0;		% Period (in seconds) of the breathing motion. Default: 7.0
        
        
        %%% Robustness simulation
        Robustness_Mode            = false; 	% Enable/Disable the robustness verification mode. Default: False
        %Scenario_selection         = 'All'		% Method for scenario selection: All (simulate all combinations), Random (randomly sample scenarios). Default: All
        Simulate_nominal_plan      = true;		% Simulate the nominal plan (without any systematic or random uncertainty). Default: True
        %Systematic_Setup_Error     = [0.25 0.25 0.25];	% Systematic error for the patient setup along the XYZ axes (cm). Default: 0.25 0.25 0.25
        %Random_Setup_Error         = [0.1  0.1  0.1];	% Standard deviation of the patient random setup error along the XYZ axes (cm). Default: 0.1 0.1 0.1
        %Systematic_Range_Error     = 3.0;		% Systematic error in percent of the proton range (%). Default: 3.0
        %Systematic_Amplitude_Error = 5.0;		% Systematic error in percent of the breathing motion amplitude for 4D simulations. Default: 5.0
        %Random_Amplitude_Error	    = 5.0;		% Random error in percent of the breathing motion amplitude for 4D simulations. Default: 5.0
        %Systematic_Period_Error	= 5.0;		% Systematic error in percent of the breathing motion period for simulations of interplay with dynamic delivery. Default: 5.0
        %Random_Period_Error        = 5.0;		% Random error in percent of the breathing motion period for simulations of interplay with dynamic delivery. Default: 5.0
        
        
        %%% Beamlet simulation
        Beamlet_Mode			= false; 	% Enable/Disable the beamlet computation mode. Default: False
        Beamlet_Parallelization = false;	% Parallelization on beamlet level is sometimes faster for beamlet simulation. This requires more memory. Default: False
        
        
        %%% Output parameters
        Output_Directory =  'MCsquareOutput';	% Name of the output directory. Default: Outputs
        
        Energy_ASCII_Output	 = false;	% Enable/Disable the output of Energy in ASCII format. Default: False
        Energy_MHD_Output    = false;	% Enable/Disable the output of Energy in MHD format. Default: False
        Energy_Sparse_Output = false;	% Enable/Disable the output of Energy in Sparse matrix format. Default: False
        Dose_ASCII_Output	 = false;	% Enable/Disable the output of Dose in ASCII format. Default: False
        Dose_MHD_Output		 = true;	% Enable/Disable the output of Dose in MHD format. Default: True
        Dose_Sparse_Output	 = true;	% Enable/Disable the output of Dose in Sparse matrix format. Default: False
        LET_ASCII_Output     = false;	% Enable/Disable the output of LET in ASCII format. Default: False
        LET_MHD_Output		 = false;	% Enable/Disable the output of LET in MHD format. Default: False
        LET_Sparse_Output	 = false;	% Enable/Disable the output of LET in Sparse matrix format. Default: False
        
        Densities_Output = false;	% Enable/Disable the export of the density map (converted from the CT image). Default: False
        Materials_Output = false;	% Enable/Disable the export of the map of materials (converted from the CT image). Default: False
        
        Compute_DVH = false;	% Enable/Disable the computation and export of DVH based on RT-Struct binary masks. Default: False
        
        Dose_Sparse_Threshold	= 0;	% The dose values above the threshold will be stored in the sparse matrix file. Default: 0
        Energy_Sparse_Threshold	= 0;	% The energy values above the threshold will be stored in the sparse matrix file. Default: 0
        LET_Sparse_Threshold	= 0;	% The LET values above the threshold will be stored in the sparse matrix file. Default: 0
        
        Score_PromptGammas	= false;	% Enable/Disable the scoring of Prompt Gammas (emitted during nuclear interactions). Default: False
        PG_LowEnergyCut     = 0.0;	% Disable the scoring of Prompt Gammas with energy below this value (MeV).  Default: 0.0
        PG_HighEnergyCut    = 50.0;	% Disable the scoring of Prompt Gammas with energy above this value (MeV).  Default: 50.0
        % Typical gamma camera would be sensitive between 3.0 and 6.0 MeV
        PG_Spectrum_NumBin  = 150;	% Number of bins to score the Prompt Gamma energy spectrum.  Default: 150
        PG_Spectrum_Binning = 0.1;	% Bin width (MeV) for the scoring of Prompt Gamma spectrum.  Default: 0.1
        
        LET_Calculation_Method	= 'StopPow'; % Select the method employed for the calculation of LET (DepositedEnergy, StopPow). Default: StopPow
        
        %Export_Beam_dose         = 'Disabled' % Export dose distribution for each beam (Enable) or entire plan (Disable). Default: Disable
        Dose_to_Water_conversion = 'Disabled'; % Select the method employed to convert simulation results (dose to medium) to dose to water (Disabled, PostProcessing, OnlineSPR). Default: Disabled
        
        Dose_Segmentation                  = false;	% Enable/Disable a segmentation of the dose map based on a density thresholding (remove dose artifacts in the air). Default: False
        Density_Threshold_for_Segmentation = 0.01;	% Density threshold employed for the segmentation (in g/cm3). Default: 0.01
    end
    
    methods
        function obj = MatRad_MCsquareConfig()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
        end
        
        function write(obj,fid)
           
            MCsquareProperties = fieldnames(obj);
            
            logicalString = {'False', 'True'};

            for i = 1:numel(MCsquareProperties)
                
                % modify fieldnames beginning with "4D"
                if strncmp(MCsquareProperties{i},'fourD',5)
                    writeString = ['4D' MCsquareProperties{i}(6:end)];
                else
                    writeString = MCsquareProperties{i};
                end

                if isa(obj.(MCsquareProperties{i}),'logical')
                    fprintf(fid,[writeString ' ' logicalString{obj.(MCsquareProperties{i})+1} '\n']);
                elseif isa(obj.(MCsquareProperties{i}),'double')
                    fprintf(fid,[writeString ' ' num2str(obj.(MCsquareProperties{i})) '\n']);
                elseif isa(obj.(MCsquareProperties{i}),'char')
                    fprintf(fid,[writeString ' ' obj.(MCsquareProperties{i}) '\n']);
                else
                    error('export not defined');
                end

            end
            
        end
        
    end
end

