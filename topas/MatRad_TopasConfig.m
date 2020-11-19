classdef MatRad_TopasConfig < handle
% MatRad_TopasConfig class definition
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
        matRad_cfg = MatRad_Config.instance(); %Instance of matRad configuration class
        
        topasExecCommand; %Defaults will be set during construction according to TOPAS installation instructions and used system
        
        parallelRuns = false; %Starts runs in parallel
                
        workingDir; %working directory for the simulation
        
        label = 'matRad_plan';
        
        %Simulation parameters
        numThreads = 0; %number of used threads, 0 = max number of threads (= num cores)
        numOfRuns = 5; %Default number of runs / batches
        modeHistories = 'num'; %'frac';
        fracHistories = 1e-4; %Fraction of histories to compute
        numHistories = 1e6; %Number of histories to compute        
        verbosity = struct( 'timefeatures',0,...
                            'cputime',true,...
                            'run',0,...
                            'event',0,...
                            'tracking',0,...
                            'material',0,...
                            'maxinterruptedhistories',1000);
        
        minRelWeight = .00001; %Threshold for discarding beamlets. 0 means all weights are being considered, can otherwise be assigned to min(w)
              
        useOrigBaseData = false; % base data of the original matRad plan will be used?
        beamProfile = 'biGaussian'; %'biGaussian' (emittance); 'simple'
        
        %Not yet implemented
        %beamletMode = false; %In beamlet mode simulation will be performed for a dose influence matrix (i.e., each beamlet simulates numHistories beamlets)
        
        pencilBeamScanning = true; %This should be always true (enables deflection)
        
        
        
        %Image
        materialConversion = 'HUToWaterSchneider'; %'CustomWaterRSP';
        %materialConversion = 'CustomWaterRSP'; 
        arrayOrdering = 'F'; %'C';
        rsp_basematerial = 'Water';
        rsp_vecLength = 10000; %Bins for RSP values when using RSP conversion with custom image converter
        rsp_methodCube = 2; %1: TsBox with variable voxels, 2: TsImageCube with Density Bins and Custom Image converter
        
                      
        
        
        %Scoring
        addVolumeScorers = true
        scoreDose = true;
        scoreTrackCount = false;
        %scoreLET = true;
        scoreReportQuantity = 'Sum'; 
        outputType = 'binary'; %'csv'; 'binary';%             
        
        %Physics 
        electronProductionCut = 0.5; %in mm
        modules_protons     = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_GenericIon  = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay','g4h-elastic_HP','g4stopping','g4ion-QMD','g4radioactivedecay'};
        modules_gamma       = {'g4em-standard_opt4','g4h-phy_QGSP_BIC_HP','g4decay'};
        
        %Geometry / World
        worldMaterial = 'G4_AIR';
        
        %filenames
        outfilenames = struct(  'patientParam','matRad_cube.txt',...
                                'patientCube','matRad_cube.dat');
                                
        infilenames = struct(   'geometry','TOPAS_matRad_geometry.txt.in',...                               
                                'surfaceScorer','TOPAS_scorer_surfaceIC.txt.in',...
                                'doseScorer','TOPAS_scorer_dose.txt.in');
        
           
    end
    
    properties (SetAccess = private)
        thisFolder;
        
        MCparam; %Struct with parameters of last simulation to be saved to file
    end
    
    methods
        function obj = MatRad_TopasConfig()
            %MatRad_MCsquareConfig Construct configuration Class for TOPAS
            %   Default execution paths are set here
            
            obj.thisFolder = fileparts(mfilename('fullpath'));
            
            obj.workingDir = [obj.thisFolder filesep 'MCrun' filesep];
            
            %Let's set some default commands taken from topas installation
            %instructions for mac & debain/ubuntu
            if ispc %We assume topas is installed in wsl (since no windows version)
                obj.topasExecCommand = 'wsl export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            elseif ismac
                obj.topasExecCommand = 'export TOPAS_G4_DATA_DIR=/Applications/G4Data; export QT_QPA_PLATFORM_PLUGIN_PATH=/Applications/topas/Frameworks; /Applications/topas/bin/topas';
            elseif isunix
                obj.topasExecCommand = 'export TOPAS_G4_DATA_DIR=~/G4Data; ~/topas/bin/topas';
            else
                obj.topasExecCommand = '';
            end               
                 
        end
        

        
        function writeAllFiles(obj,ct,pln,stf,topasBaseData,w)
            
            %Reset MCparam structure
            obj.MCparam = struct();
            obj.MCparam.tallies = {};
            obj.MCparam.nbRuns = obj.numOfRuns;
            obj.MCparam.simLabel = obj.label;
            
                        
            if ~exist(obj.workingDir,'dir')
                mkdir(obj.workingDir);
                obj.matRad_cfg.dispInfo('Created TOPAS working directory %s\n',obj.workingDir);
            end
            
            obj.MCparam.workingDir = obj.workingDir;
            
            obj.matRad_cfg.dispInfo('Writing parameter files to %s\n',obj.workingDir);

            obj.writePatient(ct,pln);
            obj.writeStfFields(ct,stf,topasBaseData,w);                       

            obj.matRad_cfg.dispInfo('Successfully written TOPAS setup files!\n')
            
            obj.writeMCparam();
        end

    end
    
    %Private sub functions for writing (private so the state of the configuration
    %can not be corrupted)
    methods (Access = private)
        
        function writeRunHeader(obj,fID,fieldIx,runIx)
           
           fprintf(fID,'s:Sim/PlanLabel = "%s"\n',obj.label);
           fprintf(fID,'s:Sim/ScoreLabel = "score_%s_field%d_run%d"\n',obj.label,fieldIx,runIx);           
           fprintf(fID,'\n');
           
           logicalString = {'"False"', '"True"'};
           
           fprintf(fID,'i:Ma/Verbosity = %d\n',obj.verbosity.material);
           fprintf(fID,'i:Ts/TrackingVerbosity = %d\n',obj.verbosity.tracking);
           fprintf(fID,'i:Ts/EventVerbosity = %d\n',obj.verbosity.event);
           fprintf(fID,'i:Ts/RunVerbosity = %d\n',obj.verbosity.run);
           fprintf(fID,'b:Ts/ShowCPUTime = %s\n',logicalString{obj.verbosity.cputime + 1});
           fprintf(fID,'i:Tf/Verbosity = %d\n',obj.verbosity.timefeatures);
           fprintf(fID,'i:Ts/MaxInterruptedHistories = %d\n',obj.verbosity.maxinterruptedhistories);
           fprintf(fID,'Ts/NumberOfThreads = %d\n',obj.numThreads);
           fprintf(fID,'i:Ts/ShowHistoryCountAtInterval = %d\n',10^(floor(log10(1/obj.numOfRuns * obj.numHistories))-1));
           fprintf(fID,'\n');
           
           
           fprintf(fID,'s:Sim/DoseScorerOutputType = "%s"\n',obj.outputType);
           fprintf(fID,'sv:Sim/DoseScorerReport = 1 "%s"\n',obj.scoreReportQuantity);           
           fprintf(fID,'\n');
           
           %fprintf(fID,'includeFile = %s/TOPAS_Simulation_Setup.txt\n',obj.thisFolder);
           %fprintf(fID,'includeFile = %s/TOPAS_matRad_geometry.txt\n',obj.thisFolder);
           %fprintf(fID,'includeFile = %s/TOPAS_scorer_surfaceIC.txt\n',obj.thisFolder);
        end
        
        function writeFieldHeader(obj,fID,fieldIx)
           fprintf(fID,'u:Sim/HalfValue = %d\n',0.5);
           fprintf(fID,'u:Sim/SIGMA2FWHM = %d\n',2.354818);
           fprintf(fID,'u:Sim/FWHM2SIGMA = %d\n',0.424661);
           fprintf(fID,'\n');
           
           fprintf(fID,'d:Sim/ElectronProductionCut = %f mm\n',obj.electronProductionCut);
           fprintf(fID,'s:Sim/WorldMaterial = "%s"\n',obj.worldMaterial);  
           fprintf(fID,'\n');
           
           fprintf(fID,'includeFile = %s\n',obj.outfilenames.patientParam);
           fprintf(fID,'\n');
                      
           fname = fullfile(obj.thisFolder,obj.infilenames.geometry);
           obj.matRad_cfg.dispInfo('Reading Geometry from %s\n',fname);
           world = fileread(fname);
           fprintf(fID,'%s\n\n',world);
        end
        
        function writeScorers(obj,fID)            
            
            obj.MCparam.outputType = obj.outputType;
            
            if obj.scoreDose
                fname = fullfile(obj.thisFolder,obj.infilenames.doseScorer);
                obj.matRad_cfg.dispDebug('Reading Dose Scorer from %s\n',fname);
                scorer = fileread(fname);                
                fprintf(fID,'%s\n',scorer);
                obj.MCparam.tallies{end+1} = 'physicalDose';
            end
                       
            if obj.addVolumeScorers
                fileList = dir(fullfile(obj.thisFolder,'TOPAS_scorer_volume_*.in'));
                for fileIx=1:length(fileList)
                    fname = fullfile(obj.thisFolder,fileList(fileIx).name);
                    obj.matRad_cfg.dispDebug('Reading Volume Scorer from %s\n',fname);
                    scorer = fileread(fname);
                    fprintf(fID,'%s\n',scorer);
                    tallyLabel = regexprep(fileList(fileIx).name,'TOPAS_scorer_volume_','');
                    tallyLabel = regexprep(tallyLabel,'.txt.in','');
                    obj.MCparam.tallies{end+1} = tallyLabel;
                end
            end

            if obj.scoreTrackCount
                fname = fullfile(obj.thisFolder,obj.infilenames.surfaceScorer);
                obj.matRad_cfg.dispDebug('Reading surface scorer from %s\n',fname);
                scorer = fileread(fname);
                fprintf(fID,'%s\n',scorer);
            end
        end
        
        function writeStfFields(obj,ct,stf,baseData,w)
            
            %snippet for future dij simulation
            %if ~exist('w','var')
            %    numBixels = sum([stf(:).totalNumOfBixels)];
            %    w = ones(numBixels,1);
            %end
            
            %Bookkeeping
            obj.MCparam.nbFields = length(stf);
            
            %Sanity check
            if numel(w) ~= sum([stf(:).totalNumOfBixels])
                obj.matRad_cfg.dispError('Given number of weights (#%d) doesn''t match bixel count in stf (#%d)',numel(w), sum([stf(:).totalNumOfBixels]));
            end
            
            nozzleToAxisDistance = baseData.nozzleToIso;           
            
            nParticlesTotalBixel = round(1e6*w);
            %nParticlesTotal = sum(nParticlesTotalBixel);
            maxParticlesBixel = 1e6*max(w(:));
            minParticlesBixel = round(max([obj.minRelWeight*maxParticlesBixel,1]));
            
            switch obj.modeHistories
                case 'num'
                    obj.fracHistories = obj.numHistories ./ sum(nParticlesTotalBixel);
                case 'frac'
                    obj.numHistories = sum(nParticlesTotalBixel);
                otherwise
                    obj.matRad_cfg.dispError('Invalid history setting!');
            end
                      
            nParticlesTotal = 0;
            
            for beamIx = 1:length(stf)
                
                SAD = stf(beamIx).SAD;
                
                sourceToNozzleDistance = SAD - nozzleToAxisDistance;
                
                %Selection of base data given the energies
                if obj.useOrigBaseData
                    [~,ixTmp,~] = intersect([ baseData.machine.data.energy], [stf.ray.energy]);
                    for i = 1:length(ixTmp)
                        selectedData(i) =  baseData.machine.data(ixTmp(i));
                    end
                    energies = [selectedData.energy];
                else
                    selectedData = [];
                    focusIndex = baseData.selectedFocus(baseData.energyIndex);
                    for i = 1:numel(focusIndex)
                        selectedData = [selectedData, baseData.monteCarloData(focusIndex(i), i)];
                    end
                    energies = [selectedData.NominalEnergy];
                end
                
                %Get Range Shifters in field ifpresent
                allRays = [stf(beamIx).ray];
                raShis = [allRays.rangeShifter];
                [~,ix] =  unique(cell2mat(squeeze(struct2cell(raShis))'),'rows');
                
                raShis = raShis(ix);
                ix = [raShis.ID] == 0;
                raShis = raShis(~ix);
                
                %Convert ID into readable string
                for r = 1:numel(raShis)
                    if isnumeric(raShis(r).ID)
                        raShis(r).topasID = ['RangeShifter' num2str(raShis(r).ID)];
                    else
                        raShis(r).topasID = ['RangeShifter' raShis(r).ID];
                    end
                end
                
                
                
                %get beamlet properties for each bixel in the stf and write
                %it into dataTOPAS
                currentBixel = 1;
                cutNumOfBixel = 0;
                nBeamParticlesTotal(beamIx) = 0;
                
                collectBixelIdx = [];
                
                %Loop over rays and then over spots on ray
                for rayIx = 1:stf(beamIx).numOfRays
                    for bixelIx = 1:stf(beamIx).numOfBixelsPerRay(rayIx)
                        
                        nCurrentParticles = nParticlesTotalBixel(currentBixel);
                        
                        % check whether there are (enough) particles for beam delivery
                        if (nCurrentParticles>minParticlesBixel)
                            
                            collectBixelIdx(end+1) = bixelIx;
                            cutNumOfBixel = cutNumOfBixel + 1;
                            bixelEnergy = stf(beamIx).ray(rayIx).energy(bixelIx);
                            [~,ixTmp,~] = intersect(energies, bixelEnergy);
                            
                            voxel_x = -stf(beamIx).ray(rayIx).rayPos_bev(3);
                            voxel_y = stf(beamIx).ray(rayIx).rayPos_bev(1);
                            
                            dataTOPAS(cutNumOfBixel).posX = -1.*voxel_x;
                            dataTOPAS(cutNumOfBixel).posY = voxel_y;
                            
                            dataTOPAS(cutNumOfBixel).current = uint32(obj.fracHistories*nCurrentParticles / obj.numOfRuns);
 
                            if obj.pencilBeamScanning
                                % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
                                % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
                                % note that Y' corresponds to the Y axis after the rotation of angleX around X axis                                
                                dataTOPAS(cutNumOfBixel).angleX = atan(dataTOPAS(cutNumOfBixel).posY / SAD);
                                dataTOPAS(cutNumOfBixel).angleY = atan(-dataTOPAS(cutNumOfBixel).posX ./ (SAD ./ cos(dataTOPAS(cutNumOfBixel).angleX)));
                                dataTOPAS(cutNumOfBixel).posX = (dataTOPAS(cutNumOfBixel).posX / SAD)*(SAD-nozzleToAxisDistance);
                                dataTOPAS(cutNumOfBixel).posY = (dataTOPAS(cutNumOfBixel).posY / SAD)*(SAD-nozzleToAxisDistance);
                            end
                            
                            if obj.useOrigBaseData
                                dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).energy;
                                dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).initFocus.SisFWHMAtIso(stf(beamIx).ray(rayIx).focusIx(bixelIx));
                                
                            else
                                dataTOPAS(cutNumOfBixel).energy = selectedData(ixTmp).MeanEnergy;
                                dataTOPAS(cutNumOfBixel).energySpread = selectedData(ixTmp).EnergySpread;
                                dataTOPAS(cutNumOfBixel).spotSize = selectedData(ixTmp).SpotSize1x;
                                dataTOPAS(cutNumOfBixel).divergence = selectedData(ixTmp).Divergence1x;
                                dataTOPAS(cutNumOfBixel).correlation = selectedData(ixTmp).Correlation1x;
                                dataTOPAS(cutNumOfBixel).focusFWHM = selectedData(ixTmp).FWHMatIso;
                            end
                            
                            %Add RangeShifterState
                            if ~isempty(raShis)
                                for r = 1:length(raShis)
                                    if stf(beamIx).ray(rayIx).rangeShifter(bixelIx).ID == raShis(r).ID
                                        raShiOut(r) = 0; %Range shifter is in beam path
                                    else
                                        raShiOut(r) = 1; %Range shifter is out of beam path / not used
                                    end
                                end
                                dataTOPAS(cutNumOfBixel).raShiOut = raShiOut;
                            end
                            
                            nBeamParticlesTotal(beamIx) = nBeamParticlesTotal(beamIx) + nCurrentParticles;
                            
                            
                        end
                        
                        currentBixel = currentBixel + 1;
                        
                    end
                end
                                
                nParticlesTotal = nParticlesTotal + nBeamParticlesTotal(beamIx);
                
                % discard data if the current has unphysical values
                idx = find([dataTOPAS.current] < 1);
                dataTOPAS(idx) = [];
                cutNumOfBixel = length(dataTOPAS(:));
                
                historyCount(beamIx) = uint32(obj.fracHistories * nBeamParticlesTotal(beamIx) / obj.numOfRuns);
                
                while sum([dataTOPAS.current]) ~= historyCount(beamIx)
                    % Randomly pick an index with the weigth given by the current
                    idx = 1:length(dataTOPAS);
                    % Note: as of Octave version 5.1.0, histcounts is not yet implemented
                    %       using histc instead for compatibility with MATLAB and Octave
                    %[~,~,R] = histcounts(rand(1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));
                    [~,R] = histc(rand(1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));
                    randIx = idx(R);
                    
                    if (sum([dataTOPAS(:).current]) > historyCount(beamIx))
                        if dataTOPAS(randIx).current > 1
                            dataTOPAS(randIx).current = dataTOPAS(randIx).current-1;
                        end
                    else
                        dataTOPAS(randIx).current = dataTOPAS(randIx).current+1;
                    end
                end
                
                historyCount(beamIx) = historyCount(beamIx) * obj.numOfRuns;  

                
                %sort dataTOPAS according to energy
                [~,ixSorted] = sort([dataTOPAS(:).energy]);
                dataTOPAS = dataTOPAS(ixSorted);
                
                %write TOPAS data base file
                fieldSetupFileName = sprintf('beamSetup_%s_field%d.txt',obj.label,beamIx);
                fileID = fopen(fullfile(obj.workingDir,fieldSetupFileName),'w');
                obj.writeFieldHeader(fileID,beamIx);                
                
                %Write modality specific info
                switch stf(beamIx).radiationMode
                    case 'protons'
                        fprintf(fileID,'s:Sim/ParticleName = "proton"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 1.0\n');
                        
                        particleA = 1;
                        particleZ = 1;
                        
                        modules = obj.modules_protons;
                        
                    case 'carbon'
                        fprintf(fileID,'s:Sim/ParticleName = "GenericIon(6,12)"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 12.0\n');
                        
                        particleA = 12;
                        particleZ = 16;
                        
                        modules = obj.modules_GenericIon;
                    %{    
                    case 'photons' %This modality is not yet used!
                           
                        fprintf(fileID,'s:Sim/ParticleName = "gamma"\n');
                        fprintf(fileID,'u:Sim/ParticleMass = 0\n');
                        
                        particleA = 0;
                        particleZ = 0;
                        
                        obj.modules_gamma;
                    %}
                    otherwise
                        obj.matRad_cfg.dispError('Invalid radiation mode %s!',stf.radiationMode)
                end
                
                fprintf(fileID,'d:Sim/GantryAngle = %f deg\n',stf(beamIx).gantryAngle);
                fprintf(fileID,'d:Sim/CouchAngle = %f deg\n',stf(beamIx).couchAngle);
                
                fprintf(fileID,'d:Tf/TimelineStart = 0. ms\n');
                fprintf(fileID,'d:Tf/TimelineEnd = %i ms\n', 10 * cutNumOfBixel);
                fprintf(fileID,'i:Tf/NumberOfSequentialTimes = %i\n', cutNumOfBixel);
                fprintf(fileID,'dv:Tf/Beam/Spot/Times = %i ', cutNumOfBixel);
                fprintf(fileID,num2str(linspace(10,cutNumOfBixel*10,cutNumOfBixel)));
                fprintf(fileID,' ms\n');
                %fprintf(fileID,'uv:Tf/Beam/Spot/Values = %i %s\n',cutNumOfBixel,num2str(collectBixelIdx));
                fprintf(fileID,'s:Tf/Beam/Energy/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Values = %i ', cutNumOfBixel);
                
                fprintf(fileID,num2str(particleA*[dataTOPAS.energy])); %Transform total energy with atomic number
                fprintf(fileID,' MeV\n');
                
                switch obj.beamProfile
                    case 'biGaussian'
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.energySpread]));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.spotSize]));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaXPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaXPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaXPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.divergence]));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.correlation]));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.spotSize]));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaYPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaYPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaYPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.divergence]));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.correlation]));
                        fprintf(fileID,'\n');
                    case 'simple'
                        fprintf(fileID,'s:Tf/Beam/FocusFWHM/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.focusFWHM]));
                        fprintf(fileID,' mm\n');
                end
                
                if obj.pencilBeamScanning
                    fprintf(fileID,'s:Tf/Beam/AngleX/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,num2str([dataTOPAS.angleX]));
                    fprintf(fileID,' rad\n');
                    fprintf(fileID,'s:Tf/Beam/AngleY/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,num2str([dataTOPAS.angleY]));
                    fprintf(fileID,' rad\n');
                end
                
                fprintf(fileID,'s:Tf/Beam/PosX/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.posX]));
                fprintf(fileID,' mm\n');
                fprintf(fileID,'s:Tf/Beam/PosY/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.posY]));
                fprintf(fileID,' mm\n');
                
                fprintf(fileID,'s:Tf/Beam/Current/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Current/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'iv:Tf/Beam/Current/Values = %i ', cutNumOfBixel);
                fprintf(fileID,num2str([dataTOPAS.current]));
                fprintf(fileID,'\n\n');
                
                %Range shifter in/out
                if ~isempty(raShis)
                    fprintf(fileID,'#Range Shifter States:\n');
                    for r = 1:numel(raShis)
                        fprintf(fileID,'s:Tf/Beam/%sOut/Function = "Step"\n',raShis(r).topasID);
                        fprintf(fileID,'dv:Tf/Beam/%sOut/Times = Tf/Beam/Spot/Times ms\n',raShis(r).topasID);
                        fprintf(fileID,'uv:Tf/Beam/%sOut/Values = %i ', raShis(r).topasID, cutNumOfBixel);
                        fprintf(fileID,num2str([dataTOPAS.raShiOut]));
                        fprintf(fileID,'\n\n');
                    end
                end
                
                
                % NozzleAxialDistance
                fprintf(fileID,'d:Ge/Nozzle/TransZ = -%f mm\n', nozzleToAxisDistance);
                if obj.pencilBeamScanning
                    fprintf(fileID,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotZ = 0.0 rad\n');
                end
                
                %Range Shifter Definition
                for r = 1:numel(raShis)
                    obj.writeRangeShifter(fileID,raShis(r),sourceToNozzleDistance);
                end
                                
                switch obj.beamProfile
                    case 'biGaussian'
                        TOPAS_beamSetup = fileread('TOPAS_beamSetup_biGaussian.txt.in');
                    case 'simple'
                        TOPAS_beamSetup = fileread('TOPAS_beamSetup_generic.txt.in');
                end
                
                fprintf(fileID,'%s\n',TOPAS_beamSetup);
                
                %translate patient according to beam isocenter
                fprintf(fileID,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-stf(beamIx).isoCenter(1));
                fprintf(fileID,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-stf(beamIx).isoCenter(2));
                fprintf(fileID,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-stf(beamIx).isoCenter(3));
                fprintf(fileID,'d:Ge/Patient/RotX=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotY=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotZ=0. deg\n');
                
                %load topas modules depending on the particle type
                fprintf(fileID,'# MODULES\n');
                moduleString = cellfun(@(s) sprintf('"%s"',s),modules,'UniformOutput',false);
                fprintf(fileID,'sv:Ph/Default/Modules = %d %s\n',length(modules),strjoin(moduleString,' '));
                
                fclose(fileID);
                %write run scripts for TOPAS
                for runIx = 1:obj.numOfRuns
                    runFileName = sprintf('%s_field%d_run%d.txt',obj.label,beamIx,runIx);
                    fileID = fopen(fullfile(obj.workingDir,runFileName),'w');
                    obj.writeRunHeader(fileID,beamIx,runIx);
                    fprintf(fileID,['i:Ts/Seed = ',num2str(runIx),'\n']);
                    fprintf(fileID,'includeFile = ./%s\n',fieldSetupFileName);
                    obj.writeScorers(fileID);
                    fclose(fileID);
                end                                
            end
            
            %Bookkeeping
            obj.MCparam.nbParticlesTotal = nParticlesTotal;
            obj.MCparam.nbHistoriesTotal = sum(historyCount);
            obj.MCparam.nbParticlesField = nBeamParticlesTotal;
            obj.MCparam.nbHistoriesField = historyCount;     
        end
        

        function writePatient(obj,ct,pln)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % matRad export CT RSP data for TOPAS simulation
            %
            % call
            %   matRad_exportCtTOPAS(ct, path, material)
            %
            % input
            %   ct:             ct cube
            %   runsPath:       path where to save the files for MC simulation
            %   basematerial:   base material to be scaled to corresponding RSP
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
           
            
            % the image cube contains the indexing of materials
            % since at the moment TOPAS does not support ushort
            % the materials should have indexes between 0 and 32767
            % therefore, the maximum length of the vector is 32768
            %answer = ''
            %while isempty(answer)
            %  prompt = 'Length of the imaging vector [2-32768]:';
            %  answer = input(prompt)
            %  if answer > 1 && answer <= 32768
            %    vecRSPlength=answer
            %  else
            %    answer = ''
            %  end
            %end
            
            medium = obj.rsp_basematerial;
            vecRSPlength=obj.rsp_vecLength;
            
            if isequal(obj.arrayOrdering,'C')
                disp('Exporting cube in C ordering...')
                permutation = [3 1 2];
            else
                disp('Exporting cube in FORTRAN ordering...')
                permutation = [2 1 3];
            end
            
                        
            cubeExport = obj.materialConversion; %'RSP'; %'HUSchneiderToWater';
            checkMaterial = false;
            rspCubeMethod = obj.rsp_methodCube; 
            
            
            paramFile = obj.outfilenames.patientParam;
            dataFile = obj.outfilenames.patientCube;
            
            %Bookkeeping
            obj.MCparam.imageCubeOrdering = obj.arrayOrdering;
            obj.MCparam.imageCubeConversionType = obj.materialConversion;
            obj.MCparam.imageCubeFile = obj.outfilenames.patientCube;
            obj.MCparam.imageCubeDim = ct.cubeDim;
            obj.MCparam.imageVoxelDimension = ct.resolution;
            
            
            nbVoxels = prod(ct.cubeDim);
            
            outfile = fullfile(obj.workingDir, paramFile);
            obj.matRad_cfg.dispInfo('Writing data to %s\n',outfile)
            fID = fopen(outfile,'w+');
           
             
            
            switch cubeExport
                case 'CustomWaterRSP'
                    rspCube = ct.cube{1};
                    rspCube = permute(rspCube,permutation); %  X,Y,Z ordering
                    
                    fbase = fopen(['materials/' medium '.txt'],'r');
                    while ~feof(fbase)
                        strLine = fgets(fbase); %# read line by line
                        fprintf(fID,'%s',strLine);
                    end
                    fclose(fbase);
                    
                    minRSP = min(rspCube(:));
                    maxRSP = max(rspCube(:));
                    
                    % to avoid zero density
                    if minRSP<1.e-6
                        minRSP=1.e-6;
                    end
                    
                    % case of homogenous water medium (i.e., RSP=1)
                    if (length(unique(ct.cube{1})) == 1) && (unique(ct.cube{1} == 1))
                        fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                        fprintf(fID,'s:Ge/Patient/Type= "TsBox"\n');
                        fprintf(fID,'s:Ge/Patient/Material = "%s"\n',medium);
                        fprintf(fID,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
                        fprintf(fID,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
                        fprintf(fID,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
                        fprintf(fID,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
                        fprintf(fID,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
                        fprintf(fID,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
                        
                        cube = NaN;
                        % otherwise
                    else
                        
                        % to avoid issues with homogenous media
                        if minRSP+1.e-6>maxRSP
                            warning('Use only one RSP value')
                            vecRSPlength = 2;
                            minRSP = 0.5*maxRSP;
                        end
                        
                        dRSP = (maxRSP-minRSP)/(vecRSPlength-1);
                        upperRSP = maxRSP+dRSP;
                        
                        ixRSP = round((rspCube-minRSP)/dRSP)+1;
                        
                        fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                        fprintf(fID,'i:Ma/%s/VariableDensityBins = %d\n',medium,vecRSPlength);
                        fprintf(fID,'u:Ma/%s/VariableDensityMin = %f\n',medium,minRSP);
                        fprintf(fID,'u:Ma/%s/VariableDensityMax = %f\n',medium,upperRSP);
                        
                        if rspCubeMethod == 1
                            fprintf(fID,'s:Ge/Patient/Type= "TsBox"\n');
                            fprintf(fID,'s:Ge/Patient/Material = "%s"\n',medium);
                            fprintf(fID,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
                            fprintf(fID,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
                            fprintf(fID,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
                            fprintf(fID,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
                            fprintf(fID,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
                            fprintf(fID,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
                            fprintf(fID,'sv:Ge/Patient/VoxelMaterials = %d\n',nbVoxels);
                            
                            voxelString = num2str(ixRSP(:)'-1,['"' medium '_VariableDensityBin_%d"\n']);
                            
                            %for ix=1:nbVoxels
                            %    fprintf(h,'"%s"\n',[ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)]);
                            %end
                            fprintf(fID,voxelString);
                            
                            if checkMaterial
                                for ix=1:nbVoxels
                                    rspMaterial{ix} = [ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)];
                                end
                                materialsUsed = unique(rspMaterial);
                                
                                %      fprintf(h,'sv:Sc/ExtractData/Material = %d ',length(materialsUsed))
                                %      for ix=1:length(materialsUsed)
                                %      fprintf(h,'"%s" ',materialsUsed{ix});
                                %      end
                                %      fprintf(h,'\n')
                            end
                            fclose(fID);
                            cube = rspCube;
                            
                        elseif rspCubeMethod == 2
                            fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                            fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                            fprintf(fID,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
                            fprintf(fID,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
                            fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                            fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                            fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "matrad"\n');
                            fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                            fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                            fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                            fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                            fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                            fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                            fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                            fclose(fID);
                            
                            % write data
                            fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                            fwrite(fID,ixRSP(:)-1,'short');
                            fclose(fID);
                            cube = rspCube;                           
                            
                        end
                    end
                    
                case 'HUToWaterSchneider'
                    huCube = int32(permute(ct.cubeHU{1},permutation));
                    
                    rspHlut = matRad_readHLUT('matRad_default.hlut');
                    densityCorrection = [];
                    for i = 1:size(rspHlut,1)-1
                        startVal = rspHlut(i,1);
                        endVal = rspHlut(i+1,1);
                        range = startVal:1:endVal-1;
                        densityCorrection(end+1:end+numel(range)) = matRad_interp1(rspHlut(:,1),rspHlut(:,2),range);
                    end
                    densityCorrection(end+1) = rspHlut(end,2); %add last missing value
                    
                    
                    %Write the Schneider Converter
                    fprintf(fID,'dv:Ge/Patient/DensityCorrection = %d %s %s\n',numel(densityCorrection),num2str(densityCorrection,'%f '),'g/cm3');
                    fprintf(fID,'iv:Ge/Patient/SchneiderHounsfieldUnitSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityOffset = 1 1\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactor = 1 0\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderDensityFactorOffset = 1 %d\n\n',-rspHlut(1,1));
                    fprintf(fID,'i:Ge/Patient/MinImagingValue = %d\n',rspHlut(1,1));
                    
                    %fprintf(h,'iv:Ge/Patient/SchneiderHUToMaterialSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
                    %fprintf(h,'sv:Ge/Patient/SchneiderElements = 2 "Hydrogen" "Oxygen"\n');
                    %fprintf(h,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 2 0.111894 0.888106\n');
                    
                    %At least include air?
                    fprintf(fID,'iv:Ge/Patient/SchneiderHUToMaterialSections = 3 %d %d %d\n',rspHlut(1,1),rspHlut(2,1),rspHlut(end,1)+1);
                    fprintf(fID,'sv:Ge/Patient/SchneiderElements = 4 "Hydrogen" "Oxygen" "Nitrogen" "Carbon"\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 4 0.0 0.23479269 0.76508170 0.00012561\n');
                    fprintf(fID,'uv:Ge/Patient/SchneiderMaterialsWeight2 = 4 0.111894 0.888106 0.0 0.0\n');
                    fprintf(fID,'dv:Ge/Patient/SchneiderMaterialMeanExcitationEnergy = 2 85.7 78.0 eV\n');
                    
                    
                    %Write the Patient
                    fprintf(fID,'s:Ge/Patient/Parent="World"\n');
                    fprintf(fID,'s:Ge/Patient/Type = "TsImageCube"\n');
                    fprintf(fID,'b:Ge/Patient/DumpImagingValues = "True"\n');
                    %fprintf(h,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
                    %fprintf(h,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
                    fprintf(fID,'s:Ge/Patient/InputDirectory = "./"\n');
                    fprintf(fID,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                    fprintf(fID,'s:Ge/Patient/ImagingtoMaterialConverter = "Schneider"\n');
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                    fprintf(fID,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                    fprintf(fID,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                    fprintf(fID,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                    fprintf(fID,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                    fprintf(fID,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                    fprintf(fID,'s:Ge/Patient/DataType  = "SHORT"\n');
                    %fprintf(h,'includeFile = HUtoMaterialSchneiderWater.txt');
                    fclose(fID);
                    
                    % write data
                    fID = fopen(fullfile(obj.workingDir, dataFile),'w');
                    fwrite(fID,huCube,'short');
                    fclose(fID);
                    cube = huCube;
                otherwise
                    obj.matRad_cfg.dispError('Material Conversion rule "%s" not implemented (yet)!\n',cubeExport);                    
            end
            
            obj.MCparam.imageCube = cube;

                           
        end
        
        function writeRangeShifter(obj,fID,rangeShifter,sourceToNozzleDistance)
            
            %Hardcoded PMMA range shifter for now
            pmma_rsp = 1.165;
            rsWidth = rangeShifter.eqThickness / pmma_rsp;
            
            fprintf(fID,'s:Ge/%s/Parent   = "Nozzle"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Type     = "TsBox"\n',rangeShifter.topasID);
            fprintf(fID,'s:Ge/%s/Material = "Lucite"\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLX      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLY      = 250 mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/HLZ      = %f  mm\n',rangeShifter.topasID,rsWidth/2);
            fprintf(fID,'d:Ge/%s/TransX   = 500 mm * Tf/Beam/%sOut/Value\n',rangeShifter.topasID,rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransY   = 0   mm\n',rangeShifter.topasID);
            fprintf(fID,'d:Ge/%s/TransZ   = %f mm\n',rangeShifter.topasID,rangeShifter.sourceRashiDistance - sourceToNozzleDistance);
            
        end
        
        function writeMCparam(obj)
            %write MCparam file with basic parameters
              
            MCparam = obj.MCparam;            
            save(fullfile(obj.workingDir,'MCparam.mat'),'MCparam','-v7');
        end
        
    end
end

