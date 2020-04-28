classdef MatRad_TopasBaseData < MatRad_MCemittanceBaseData
    % MatRad_TopasBaseData class for calculating TOPAS base data and
    % writing it into a file, formatted for TOPAS to use
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
        TopasConfig     %TOPAS specific data defining simulation parameters
    end
    
    methods
        function obj = MatRad_TopasBaseData(varargin)
            % set default values for the simulation, if not specified
            if  nargin < 2 % input MatRad_TopasBaseData(machine)
                TopasConfig = struct;  
                constArguments = 1;
            elseif nargin < 3 && isfield(varargin{2}, 'gantryAngle') % input MatRad_TopasBaseData(machine,stf)
                TopasConfig = struct;
                constArguments = [1,2];
            else % input MatRad_TopasBaseData(machine,T) || (machine,stf,T)
                TopasConfig = varargin{end};
                if nargin == 2
                    constArguments = 1;
                elseif nargin == 3
                    constArguments = [1,2];
                end
            end
            
            % Call MatRad_MCemmitanceBaseData constructor
            obj = obj@MatRad_MCemittanceBaseData(varargin{constArguments});
            
            obj.TopasConfig = TopasConfig;
            
            if ~isfield(obj.TopasConfig,'filepath')
                % MatRad filepath
                % default: 'topas/MCexport/'
                obj.TopasConfig.filepath = 'topas/MCexport/';
            end
            
            if ~isfield(obj.TopasConfig,'fracHistories')
                % default: 1e-4
                obj.TopasConfig.fracHistories = 1e-4;
            end
            
            if ~isfield(obj.TopasConfig,'numOfThreads')
                % default: 0, all available
                obj.TopasConfig.numOfThreads = 0;
            end
            
            if ~isfield(obj.TopasConfig,'outputType')
                % default: csv
                obj.TopasConfig.outputType = 'csv';
            end
            
            if ~isfield(obj.TopasConfig,'minRelWeight')
                % default: 0.001
                % minRelWeight = 0 means all weights are being considered
                % can otherwise be assigned to min(w)
                obj.TopasConfig.minRelWeight = 0.001;
            end
            
            if ~isfield(obj.TopasConfig,'numOfRuns')
                % default: 5 runs
                obj.TopasConfig.numOfRuns = 5;
            end
            
            if ~isfield(obj.TopasConfig,'useOrigBaseData')
                % base data of the original matRad plan will be used
                % default: true
                obj.TopasConfig.useOrigBaseData = false;
            end
            
            if ~isfield(obj.TopasConfig,'beamProfile')
                % default: 5 runs
                obj.TopasConfig.beamProfile = 'biGaussian';
            end
            
            if ~isfield(obj.TopasConfig,'pencilBeamScanning')
                %   beamProfile: 'simple' of 'biGaussian'
                % default: true
                obj.TopasConfig.pencilBeamScanning = 'true';
            end
            
            if ~isfield(obj.TopasConfig,'electronProdCut')
                % default: 0.5
                obj.TopasConfig.electronProdCut = '0.5';
            end
                      
            % Check if contradicting settings are being used
            if obj.TopasConfig.useOrigBaseData == true && ~strcmp(obj.TopasConfig.beamProfile,'simple')
                error('Original base data only usable with simple beam geometry!')
            end
        end
        
        function obj = writeTopasData(obj,ct,stf,pln,w)
            %function that writes a data file containing stf specific data
            %for a Monte Carlo simulation with TOPAS
            
            % load default modules
            if ~isfield(obj.TopasConfig,'modules')
                if strcmp(pln.radiationMode,'protons')
                    obj.TopasConfig.modules = {'"g4em-standard_opt4"','"g4h-phy_QGSP_BIC_HP"','"g4decay"','"g4ion-QMD"','"g4h-elastic_HP"','"g4stopping"','"g4radioactivedecay"'};
                else
                    obj.TopasConfig.modules = {'"g4em-standard_opt4"','"g4h-phy_QGSP_BIC_HP"','"g4decay"','"g4h-elastic_HP"','"g4stopping"','"g4ion-inclxx"'};
                end
            end
            
            %look up focus indices
            focusIndex = obj.selectedFocus(obj.energyIndex);
            
            % NozzleAxialDistance
            nozzleAxialDistance_mm = 1500;
            SAD_mm  = obj.machine.meta.SAD;
            if isfield( obj.machine.meta,'nozzleAxialDistance')
                disp('Using NAD from basedata')
                nozzleAxialDistance_mm =  obj.machine.meta.nozzleAxialDistance;
            else
                disp('Using default nozzleAxialDistance')
            end
            
            for beamIx = 1:length(stf)
                
                if obj.TopasConfig.useOrigBaseData
                    [~,ixTmp,~] = intersect([ obj.machine.data.energy], [stf.ray.energy]);
                    for i = 1:length(ixTmp)
                        selectedData(i) =  obj.machine.data(ixTmp(i));
                    end
                    energies = [selectedData.energy];
                else
                    selectedData = [];
                    for i = 1:numel(focusIndex)
                        selectedData = [selectedData, obj.monteCarloData(focusIndex(i), i)];
                    end
                    energies = [selectedData.NominalEnergy];
                end
                
                %get beamlet properties for each bixel in the stf and write
                %it into dataTOPAS
                currentBixel = 1;
                cutNumOfBixel = 0;
                nbParticlesTotal = 0;
                
                for rayIx = 1:stf(beamIx).numOfRays
                    
                    for bixelIx = 1:stf(beamIx).numOfBixelsPerRay(rayIx)
                        
                        voxel_nbParticles = round(1e6*w(currentBixel));
                        maxParticlesInSpot = 1e6*max(w(:));
                        minNbParticlesSpot = round(max([obj.TopasConfig.minRelWeight*maxParticlesInSpot,1]));
                        
                        % check whether there are (enough) particles for beam delivery
                        if (voxel_nbParticles>minNbParticlesSpot)
                            
                            cutNumOfBixel = cutNumOfBixel + 1;
                            bixelEnergy = stf(beamIx).ray(rayIx).energy(bixelIx);
                            [~,ixTmp,~] = intersect(energies, bixelEnergy);
                            
                            
                            voxel_x = -stf(beamIx).ray(rayIx).rayPos_bev(3);
                            voxel_y = stf(beamIx).ray(rayIx).rayPos_bev(1);
                            
                            dataTOPAS(cutNumOfBixel).posX = -1.*voxel_x;
                            dataTOPAS(cutNumOfBixel).posY = voxel_y;
                            
                            dataTOPAS(cutNumOfBixel).current = uint32(obj.TopasConfig.fracHistories*voxel_nbParticles);
                            
                            if obj.TopasConfig.pencilBeamScanning
                                % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
                                % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
                                % note that Y' corresponds to the Y axis after the rotation of angleX around X axis
                                dataTOPAS(cutNumOfBixel).angleX = atan(dataTOPAS(cutNumOfBixel).posY / SAD_mm);
                                dataTOPAS(cutNumOfBixel).angleY = atan(-dataTOPAS(cutNumOfBixel).posX ./ (SAD_mm ./ cos(dataTOPAS(cutNumOfBixel).angleX)));
                                dataTOPAS(cutNumOfBixel).posX = (dataTOPAS(cutNumOfBixel).posX / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);
                                dataTOPAS(cutNumOfBixel).posY = (dataTOPAS(cutNumOfBixel).posY / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);
                            end
                            
                            if obj.TopasConfig.useOrigBaseData
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
                            nbParticlesTotal = nbParticlesTotal + voxel_nbParticles;
                        end
                        
                        currentBixel = currentBixel + 1;
                        
                    end
                end
                
                % discard data if the current has unphysical values
                idx = find([dataTOPAS.current] < 1);
                dataTOPAS(idx) = [];
                
                historyCount = uint32(obj.TopasConfig.fracHistories * nbParticlesTotal);

                while sum([dataTOPAS.current]) ~= historyCount
                    % Randomly pick an index with the weigth given by the current
                    idx = 1:length(dataTOPAS);
                    [~,~,R] = histcounts(rand(1),cumsum([0;double(transpose([dataTOPAS(:).current]))./double(sum([dataTOPAS(:).current]))]));
                    randIx = idx(R);
                    
                    if (sum([dataTOPAS(:).current]) > historyCount)
                        if dataTOPAS(randIx).current > 1
                            dataTOPAS(randIx).current = dataTOPAS(randIx).current-1;
                        end
                    else
                        dataTOPAS(randIx).current = dataTOPAS(randIx).current+1;
                    end
                end
                
                
                %sort dataTOPAS according to energy
                [~,ixSorted] = sort([dataTOPAS(:).energy]);
                dataTOPAS = dataTOPAS(ixSorted);
                
                %write TOPAS data base file
                fileID = fopen([obj.TopasConfig.filepath,'beamSetup_matRad_plan_field',num2str(beamIx),'.txt'],'w');
                
                fprintf(fileID,'i:Ts/ShowHistoryCountAtInterval = %i\n',historyCount/20);
                fprintf(fileID,'s:Sim/PlanLabel = "simData_matRad_plan_field1_run" + Ts/Seed\n');
                fprintf(fileID,'d:Sim/GantryAngle = %.6f deg\n', stf(beamIx).gantryAngle);
                fprintf(fileID,'d:Sim/CouchAngle = %.6f deg\n', stf(beamIx).couchAngle);
                fprintf(fileID,'s:Sim/ParticleName = "proton"\n');
                fprintf(fileID,'u:Sim/ParticleMass = 1.0\n');
                fprintf(fileID,'i:Sim/NbThreads = %i\n', obj.TopasConfig.numOfThreads);
                fprintf(fileID,'d:Tf/TimelineStart = 0. ms\n');
                fprintf(fileID,'d:Tf/TimelineEnd = %i ms\n', 10 * cutNumOfBixel);
                fprintf(fileID,'i:Tf/NumberOfSequentialTimes = %i\n', cutNumOfBixel);
                fprintf(fileID,'dv:Tf/Beam/Spot/Times = %i ', cutNumOfBixel);
                fprintf(fileID,strjoin(string(linspace(10,cutNumOfBixel*10,cutNumOfBixel))));
                fprintf(fileID,' ms\n');
                fprintf(fileID,'s:Tf/Beam/Energy/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/Energy/Values = %i ', cutNumOfBixel);
                fprintf(fileID,strjoin(string([dataTOPAS.energy])));
                fprintf(fileID,' MeV\n');
                
                switch obj.TopasConfig.beamProfile
                    case 'biGaussian'
                        fprintf(fileID,'s:Tf/Beam/EnergySpread/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/EnergySpread/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/EnergySpread/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.energySpread])));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.spotSize])));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaXPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaXPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaXPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.divergence])));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationX/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationX/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationX/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.correlation])));
                        fprintf(fileID,'\n');
                        
                        fprintf(fileID,'s:Tf/Beam/SigmaY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.spotSize])));
                        fprintf(fileID,' mm\n');
                        fprintf(fileID,'s:Tf/Beam/SigmaYPrime/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/SigmaYPrime/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/SigmaYPrime/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.divergence])));
                        fprintf(fileID,'\n');
                        fprintf(fileID,'s:Tf/Beam/CorrelationY/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/CorrelationY/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'uv:Tf/Beam/CorrelationY/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.correlation])));
                        fprintf(fileID,'\n');
                    case 'simple'
                        fprintf(fileID,'s:Tf/Beam/FocusFWHM/Function = "Step"\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Times = Tf/Beam/Spot/Times ms\n');
                        fprintf(fileID,'dv:Tf/Beam/FocusFWHM/Values = %i ', cutNumOfBixel);
                        fprintf(fileID,strjoin(string([dataTOPAS.focusFWHM])));
                        fprintf(fileID,' mm\n');
                end
                
                if obj.TopasConfig.pencilBeamScanning
                    fprintf(fileID,'s:Tf/Beam/AngleX/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleX/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,strjoin(string([dataTOPAS.angleX])));
                    fprintf(fileID,' rad\n');
                    fprintf(fileID,'s:Tf/Beam/AngleY/Function = "Step"\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Times = Tf/Beam/Spot/Times ms\n');
                    fprintf(fileID,'dv:Tf/Beam/AngleY/Values = %i ', cutNumOfBixel);
                    fprintf(fileID,strjoin(string([dataTOPAS.angleY])));
                    fprintf(fileID,' rad\n');
                end
                
                fprintf(fileID,'s:Tf/Beam/PosX/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosX/Values = %i ', cutNumOfBixel);
                fprintf(fileID,strjoin(string([dataTOPAS.posX])));
                fprintf(fileID,' mm\n');
                fprintf(fileID,'s:Tf/Beam/PosY/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'dv:Tf/Beam/PosY/Values = %i ', cutNumOfBixel);
                fprintf(fileID,strjoin(string([dataTOPAS.posY])));
                fprintf(fileID,' mm\n');
                
                fprintf(fileID,'s:Tf/Beam/Current/Function = "Step"\n');
                fprintf(fileID,'dv:Tf/Beam/Current/Times = Tf/Beam/Spot/Times ms\n');
                fprintf(fileID,'iv:Tf/Beam/Current/Values = %i ', cutNumOfBixel);
                fprintf(fileID,strjoin(string(round([dataTOPAS.current]))));
                fprintf(fileID,'\n\n');
                                
                if isfield(pln,'propStf')
                    fprintf(fileID,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-pln.propStf.isoCenter(beamIx,1));
                    fprintf(fileID,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-pln.propStf.isoCenter(beamIx,2));
                    fprintf(fileID,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-pln.propStf.isoCenter(beamIx,3));
                else
                    fprintf(fileID,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-pln.isoCenter(beamIx,1));
                    fprintf(fileID,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-pln.isoCenter(beamIx,2));
                    fprintf(fileID,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-pln.isoCenter(beamIx,3));
                end
                fprintf(fileID,'d:Ge/Patient/RotX=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotY=0. deg\n');
                fprintf(fileID,'d:Ge/Patient/RotZ=0. deg\n');
                
                fprintf(fileID,'includeFile = ./matRad_RSPcube.txt\n');
                fprintf(fileID,'###################\n');
                
                % NozzleAxialDistance
                fprintf(fileID,'d:Ge/Nozzle/TransZ = -%f mm\n', nozzleAxialDistance_mm);
                if obj.TopasConfig.pencilBeamScanning
                    fprintf(fileID,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotZ = 0.0 rad\n');
                end
                
                fprintf(fileID,['d:Ph/Default/CutForElectron = ',obj.TopasConfig.electronProdCut,' mm\n']);
                
                switch obj.TopasConfig.beamProfile
                    case 'biGaussian'
                        TOPAS_beamSetup = fopen(['TOPAS_beamSetup_biGaussian_' pln.radiationMode '.txt'],'r');
                    case 'simple'
                        TOPAS_beamSetup = fopen(['TOPAS_beamSetup_generic_' pln.radiationMode '.txt'],'r');
                end
                
                % copy standard values from TOPAS_beamSetup
                while ~feof(TOPAS_beamSetup)
                    strLine = fgets(TOPAS_beamSetup); %# read line by line
                    fprintf(fileID,'%s',strLine);
                end
                
                fprintf(fileID,'\n');
                fprintf(fileID,['s:Sc/Patient/Tally_DoseToMedium/OutputType = "', obj.TopasConfig.outputType ,'"\n']);
                fprintf(fileID,'\n');
          
                %load topas modules depending on the particle type
                fprintf(fileID,'# MODULES\n');
                fprintf(fileID,['sv:Ph/Default/Modules = ',num2str(length(obj.TopasConfig.modules)),' ', strjoin(obj.TopasConfig.modules) ,'\n']);
                
                fclose(fileID);
                %write run scripts for TOPAS
                basematerial = '';
                if ~exist('machine') || ~isfield( obj.machine.meta,'basematerial')
                    disp('Using default base material: Water')
                    basematerial = 'Water';
                else
                    basematerial =  obj.machine.meta.basematerial;
                end
                
                rspCube = matRad_exportCtTOPAS(ct, obj.TopasConfig.filepath, basematerial);
                
                for runIx = 1:obj.TopasConfig.numOfRuns
                    fileID = fopen([obj.TopasConfig.filepath,'matRad_plan_field',num2str(beamIx),'_run',num2str(runIx),'.txt'],'w');
                    fprintf(fileID,['i:Ts/Seed = ',num2str(runIx),'\n']);
                    fprintf(fileID,'includeFile = ./beamSetup_matRad_plan_field1.txt');
                    fclose(fileID);
                end
                
                %write MCparam file with basic parameters
                MCparam.nbRuns = obj.TopasConfig.numOfRuns;
                MCparam.nbFields = length(stf);
                MCparam.tallies = {'physicalDose'};
                MCparam.simLabel = 'matRad_plan';
                MCparam.cubeDim = ct.cubeDim;
                MCparam.RSP = rspCube;
                MCparam.voxelDimensions = ct.resolution;
                MCparam.outputType = obj.TopasConfig.outputType;
                MCparam.nbHistories = historyCount;
                MCparam.nbParticles = nbParticlesTotal;
                save([obj.TopasConfig.filepath,'MCparam.mat'],'MCparam');
                
                disp('Successfully written TOPAS setup files!')
            end
        end
    end
end

