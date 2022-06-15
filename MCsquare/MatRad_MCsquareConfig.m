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
        %%% Parameter for continuity
        engine = 'MCsquare';
        
        %%% Simulation parameters:
        Num_Threads   =	0;		% Number of parallel calculation threads. Default: 0 = max available threads
        RNG_Seed      =	0;		% Seed for the random number generator (deterministic result only with single thread). Default: 0 = seed based on the time
        numHistories  = 1e6;		% Number of primary protons to simulate. Default: 1e7
        E_Cut_Pro     =	0.5;		% Energy cut (in MeV) below which heavy charged particles are locally absorbed. Default: 0.5
        D_Max	      =	0.2;		% Maximum distance between two step (cm). Default: 0.2
        Epsilon_Max   =	0.25;		% Fractional energy loss (dE/T) per step. Default: 0.25
        Te_Min	      =	0.05;		% Threshold energy (MeV) for the production of secondary electrons (currently locally absorbed). Default: 0.05
        Stat_uncertainty    = 0.0	% Maximum statistical uncertainty (in percent). Default: 0.0 = no maximum uncertainty (number of proton = Num_Primaries)
        % As a reference: 200 MeV protons can transfer a maximum energy of 0.5 MeV to ?-electrons which correspond to a range of 7 mm in lung tissues.
        
        %%% Input files
        CT_File                     = 'Patient.mhd';				% Name of the CT file. Default: CT.mhd
        HU_Density_Conversion_File	= 'Scanners/matRad_default/HU_Density_Conversion.txt';	% Name of the file containing HU to density conversion data. Default: HU_Density_Conversion.txt
        HU_Material_Conversion_File	= 'Scanners/matRad_default/HU_Material_Conversion.txt';	% Name of the file containing HU to material conversion data. Default: HU_Material_Conversion.txt
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
            %MatRad_MCsquareConfig Configuration Class for MCsquare           
        end

        function writeMCsquareinputAllFiles(obj,filename,stf)
            % generate input files for MCsquare dose calcualtion from matRad
            %
            % call
            %   obj.writeMCsquareinputAllFiles(filename,MCsquareConfig,stf)
            %
            % input
            %   filename:       filename
            %   stf:            matRad steering information struct
            %
            % output
            %   -
            %
            % References
            %   [1] https://openreggui.org/git/open/REGGUI/blob/master/functions/io/convert_Plan_PBS.m
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


            %% write overall configuration file
            fileHandle = fopen(filename,'w');
            obj.write(fileHandle);
            fclose(fileHandle);

            %% prepare steering file writing
            numOfFields = length(stf);
            if obj.Beamlet_Mode
                totalMetersetWeightOfAllFields = 1;
            else
                totalMetersetWeightOfFields = NaN*ones(numOfFields,1);
                for i = 1:numOfFields
                    totalMetersetWeightOfFields(i) = sum([stf(i).energyLayer.numOfPrimaries]);
                end
                totalMetersetWeightOfAllFields = sum(totalMetersetWeightOfFields);
            end

            %% write steering file

            fileHandle = fopen(obj.BDL_Plan_File,'w');

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
            fprintf(fileHandle,[num2str(totalMetersetWeightOfAllFields) '\n']);

            for i = 1:numOfFields
                fprintf(fileHandle,'\n#FIELD-DESCRIPTION\n');
                fprintf(fileHandle,'###FieldID\n');
                fprintf(fileHandle,[num2str(i) '\n']);
                fprintf(fileHandle,'###FinalCumulativeMeterSetWeight\n');
                if obj.Beamlet_Mode
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

                %Range shfiter
                if stf(i).rangeShifterID ~= 0
                    fprintf(fileHandle,'###RangeShifterID\n%d\n',stf(i).rangeShifterID);
                    fprintf(fileHandle,'###RangeShifterType\n%s\n',stf(i).rangeShifterType);
                end

                metersetOffset = 0;
                fprintf(fileHandle,'\n#SPOTS-DESCRIPTION\n');
                for j = 1:numOfEnergies
                    fprintf(fileHandle,'####ControlPointIndex\n');
                    fprintf(fileHandle,[num2str(j) '\n']);
                    fprintf(fileHandle,'####SpotTunnedID\n');
                    fprintf(fileHandle,['1\n']);
                    fprintf(fileHandle,'####CumulativeMetersetWeight\n');
                    if obj.Beamlet_Mode
                        cumulativeMetersetWeight = j/numOfEnergies * 1/numOfFields;
                    else
                        cumulativeMetersetWeight = metersetOffset + sum([stf(i).energyLayer(j).numOfPrimaries]);
                        metersetOffset = cumulativeMetersetWeight;
                    end
                    fprintf(fileHandle,[num2str(cumulativeMetersetWeight) '\n']);
                    fprintf(fileHandle,'####Energy (MeV)\n');
                    fprintf(fileHandle,[num2str(stf(i).energies(j)) '\n']);

                    %Range shfiter
                    if stf(i).rangeShifterID ~= 0
                        rangeShifter = stf(i).energyLayer(j).rangeShifter;
                        if rangeShifter.ID ~= 0
                            fprintf(fileHandle,'####RangeShifterSetting\n%s\n','IN');
                            pmma_rsp = 1.165; %TODO: hardcoded for now
                            rsWidth = rangeShifter.eqThickness / pmma_rsp;
                            isoToRaShi = stf(i).SAD - rangeShifter.sourceRashiDistance + rsWidth;
                            fprintf(fileHandle,'####IsocenterToRangeShifterDistance\n%f\n',-isoToRaShi/10); %in cm
                            fprintf(fileHandle,'####RangeShifterWaterEquivalentThickness\n%f\n',rangeShifter.eqThickness);
                        else
                            fprintf(fileHandle,'####RangeShifterSetting\n%s\n','OUT');
                        end
                    end

                    fprintf(fileHandle,'####NbOfScannedSpots\n');
                    numOfSpots = size(stf(i).energyLayer(j).targetPoints,1);
                    fprintf(fileHandle,[num2str(numOfSpots) '\n']);
                    fprintf(fileHandle,'####X Y Weight\n');
                    for k = 1:numOfSpots
                        if obj.Beamlet_Mode
                            n = stf(i).energyLayer(j).numOfPrimaries(k);
                        else
                            n = stf(i).energyLayer(j).numOfPrimaries(k) / obj.mcSquare_magicFudge(stf(i).energies(j));
                        end
                        fprintf(fileHandle,[num2str(stf(i).energyLayer(j).targetPoints(k,:)) ' ' num2str(n) '\n']);
                    end
                end
            end

            fclose(fileHandle);

        end

        function gain = mcSquare_magicFudge(~,energy)
            % mcSquare will scale the spot intensities in
            % https://gitlab.com/openmcsquare/MCsquare/blob/master/src/data_beam_model.c#L906
            % by this factor so we need to divide up front to make things work. The
            % original code can be found at https://gitlab.com/openmcsquare/MCsquare/blob/master/src/compute_beam_model.c#L16

            K = 35.87; % in eV (other value 34.23 ?)

            % // Air stopping power (fit ICRU) multiplied by air density
            SP = (9.6139e-9*energy^4 - 7.0508e-6*energy^3 + 2.0028e-3*energy^2 - 2.7615e-1*energy + 2.0082e1) * 1.20479E-3 * 1E6; % // in eV / cm

            % // Temp & Pressure correction
            PTP = 1.0;

            % // MU calibration (1 MU = 3 nC/cm)
            % // 1cm de gap effectif
            C = 3.0E-9; % // in C / cm

            % // Gain: 1eV = 1.602176E-19 J
            gain = (C*K) / (SP*PTP*1.602176E-19);

            % divide by 1e7 to not get tiny numbers...
            gain = gain/1e7;

        end

        function writeMhd(obj,cube,resolution)
            % References
            %   -
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2020 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% write header file
            fileHandle = fopen(obj.CT_File,'w');

            fprintf(fileHandle,'ObjectType = Image\n');
            fprintf(fileHandle,'NDims = 3\n');
            fprintf(fileHandle,'BinaryData = True\n');
            fprintf(fileHandle,'BinaryDataByteOrderMSB = False\n');
            fprintf(fileHandle,'CompressedData = False\n');
            fprintf(fileHandle,'TransformMatrix = 1 0 0 0 1 0 0 0 1\n');
            fprintf(fileHandle,'Offset = 0 0 0\n');
            fprintf(fileHandle,'CenterOfRotation = 0 0 0\n');
            fprintf(fileHandle,'AnatomicalOrientation = RAI\n');
            fprintf(fileHandle,'ElementSpacing = %f %f %f\n',resolution);
            fprintf(fileHandle,'DimSize = %d %d %d\n',size(cube,2),size(cube,1),size(cube,3));
            fprintf(fileHandle,'ElementType = MET_DOUBLE\n');
            filenameRaw = [obj.CT_File(1:end-4) '.raw'];
            fprintf(fileHandle,'ElementDataFile = %s\n',filenameRaw);

            fclose(fileHandle);

            %% write data file
            dataFileHandle = fopen(filenameRaw,'w');

            cube = flip(cube,2);
            cube = permute(cube,[2 1 3]);

            fwrite(dataFileHandle,cube(:),'double');
            fclose(dataFileHandle);
        end

        function cube = readMhd(~,filename)
            % matRad mhd file reader
            %
            % call
            %   cube = matRad_readMhd(folder,filename)
            %
            % input
            %   folder:   folder where the *raw and *mhd file are located
            %   filename: filename
            %
            % output
            %   cube:     3D array
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


            %% read header
            headerFileHandle = fopen([obj.Output_Directory, filesep filename],'r');

            s = textscan(headerFileHandle, '%s', 'delimiter', '\n');

            % read dimensions
            idx = find(~cellfun(@isempty,strfind(s{1}, 'DimSize')),1,'first');
            dimensions = cell2mat(textscan(s{1}{idx},'DimSize = %f %f %f'));

            % read filename of data
            idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementDataFile')),1,'first');
            tmp = textscan(s{1}{idx},'ElementDataFile = %s');
            dataFilename = cell2mat(tmp{1});

            % get data type
            idx = find(~cellfun(@isempty,strfind(s{1}, 'ElementType')),1,'first');
            tmp = textscan(s{1}{idx},'ElementType = MET_%s');
            type = lower(cell2mat(tmp{1}));

            fclose(headerFileHandle);

            %% read data
            dataFileHandle = fopen([obj.Output_Directory filesep dataFilename],'r');
            cube = reshape(fread(dataFileHandle,inf,type),dimensions);
            cube = permute(cube,[2 1 3]);
            cube = flip(cube,2);
            fclose(dataFileHandle);
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

