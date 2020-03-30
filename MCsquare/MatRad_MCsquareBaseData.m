classdef MatRad_MCsquareBaseData
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matRad_MCsquareBaseData calculates MonteCarlo base data and creates
    % base data file either formatted for use by MCsquare or Topas
    %
    %
    %
    %
    % References
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
    
    properties
        machine         %matRad base data machine struct
        bdl_path = ''   %stores path to generated file
        nozzleToIso     %Nozzle to Isocenter Distance
        smx             %Scanning magnet X to isocenter Distance
        smy             %Scanning magnet y to isocenter Distance
        mcSquareData    %MCsquare Phase space data struct
        selectedFocus   %array containing selected focus indices per energy
        FWHMatIso       %array containing FWHM values at iscenter for every energy
    end
    
    properties (SetAccess = private)
        stfCompressed   %measure whether function has additional info about
        %the stf
        problemSigma    % = 1, when there was a problem calculating sigma
        energyIndex     %Indices of calculated energies
    end
    
    methods
        function obj = MatRad_MCsquareBaseData(machine,stf)
            %MatRad_MCsquareBaseData construct an instance of the MCsquare
            %Base data format using a focus index
            
            %stfCompressed states whether mcSquareData are calculated for
            %all energies (false) or only for energies which exist in given
            %stf. If function is called without stf stfCompressed = false.
            if nargin < 2
                obj.stfCompressed = false;
            else
                obj.stfCompressed = true;
            end
            
            obj.machine = machine;
            obj.problemSigma = false;
            obj.selectedFocus = ones(numel(machine.data),1) * NaN;
            
            if isfield(machine.meta,'BAMStoIsoDist')
                obj.nozzleToIso = machine.meta.BAMStoIsoDist;
            else
                warning('No information on BAMS to isocenter distance. Using generic value of 500mm');
                obj.nozzleToIso = 500;
            end
            
            SAD = machine.meta.SAD;
            
            obj.smx = SAD;
            obj.smy = SAD;
            
            obj.mcSquareData = [];
            
            %select needed energies and according focus indices by using stf
            if obj.stfCompressed
                tmp = [stf(:).ray];
                plannedEnergies     = [tmp.energy];
                focusIndex          = [tmp.focusIx];
                [~, ind]            = unique(plannedEnergies);
                plannedEnergies     = plannedEnergies(ind);
                focusIndex          = focusIndex(ind);
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
                
                %if no stf was refered all energies are chosen, while setting
                %the focus index for all energies to preliminary 1
            else
                plannedEnergies = [machine.data(:).energy];
                focusIndex = ones(size(plannedEnergies));
                [~ ,obj.energyIndex, ~] = intersect([machine.data(:).energy],plannedEnergies);
            end
            
            obj.selectedFocus(obj.energyIndex) = focusIndex;
            
            count = 1;
            for ii = 1:numel(obj.energyIndex)
                
                i = obj.energyIndex(ii);
                
                %look up whether MonteCarlo data are already present in
                %machine file , if so do not recalculate
                if isfield(machine.data(i),'mcSquareData')
                    if (isempty(machine.data(i).mcSquareData) == 0)
                        obj.mcSquareData = [obj.mcSquareData, machine.data(i).mcSquareData];
                        count = count + 1;
                        continue;
                    end
                end
                
                
                %calculate mcSquareData for given energy and every focus
                %index
                data = [];
                energyData = obj.fitPhaseSpaceForEnergy(i);
                obj.FWHMatIso = [];
                for j = 1:size(machine.data(i).initFocus.sigma,1)
                    
                    tmp = energyData;
                    opticsData = obj.fitBeamOpticsForEnergy(i, j);
                    
                    f = fieldnames(opticsData);
                    for a = 1:length(f)
                        tmp.(f{a}) = opticsData.(f{a});
                    end
                    
                    data = [data; tmp];
                end
                
                obj.mcSquareData = [obj.mcSquareData, data];
                
                count = count + 1;
            end
            
            %throw out warning if there was a problem in calculating the
            %width of the Bragg peak in obj.fitBeamOpticsForEnergy
            if obj.problemSigma
                warning('Calculation of FWHM of bragg peak in base data not possible! Using simple approximation for energy spread');
            end
        end
        
        function mcDataEnergy = fitPhaseSpaceForEnergy(obj,energyIx)
            %function to calculate mean energy and energy spread used by
            %mcSquare for given energy
            
            i = energyIx;
            
            mcDataEnergy.NominalEnergy = obj.machine.data(i).energy;
            
            newDepths = linspace(0,obj.machine.data(i).depths(end),numel(obj.machine.data(i).depths) * 100);
            newDose   = interp1(obj.machine.data(i).depths, obj.machine.data(i).Z, newDepths, 'spline');
            
            %find FWHM w50 of bragg peak and range of 80% does fall off
            [maxV, maxI] = max(newDose);
            [~, r80ind] = min(abs(newDose(maxI:end) - 0.8 * maxV));
            r80ind = r80ind - 1;
            r80 = interp1(newDose(maxI + r80ind - 1:maxI + r80ind + 1), ...
                newDepths(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV) ...
                + obj.machine.data(i).offset;
            
            
            [~, d50rInd] = min(abs(newDose(maxI:end) - 0.5 * maxV));
            d50rInd = d50rInd - 1;
            d50_r = interp1(newDose(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                newDepths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
            
            if (newDose(1) < 0.5 * maxV)
                [~, d50lInd] = min(abs(newDose(1:maxI) - 0.5*maxV));
                d50_l = interp1(newDose(d50lInd - 1:d50lInd + 1), ...
                    newDepths(d50lInd - 1:d50lInd + 1), 0.5 * maxV);
                w50 = d50_r - d50_l;
                %if width left of peak cannot be determined use r80 as width
            else
                d50_l = newDepths(maxI);
                w50 = r80;
                obj.problemSigma = true;
            end
            
            %calcualte mean energy used my mcSquare with a formula fitted
            %to TOPAS data
            meanEnergy = @(x) 5.762374661332111e-20 * x^9 - 9.645413625310569e-17 * x^8 + 7.073049219034644e-14 * x^7 ...
                - 2.992344292008054e-11 * x^6 + 8.104111934547256e-09 * x^5 - 1.477860913846939e-06 * x^4 ...
                + 1.873625800704108e-04 * x^3 - 1.739424343114980e-02 * x^2 + 1.743224692623838e+00 * x ...
                + 1.827112816899668e+01;
            mcDataEnergy.MeanEnergy = meanEnergy(r80);
            
            %calculate energy straggling using formulae deducted from paper
            %"An analytical approximation of the Bragg curve for therapeutic
            %proton beams" by T. Bortfeld et al.
            totalSigmaSq = ((w50) / 6.14)^2;
            
            totalSpreadSq = @(x) 2.713311945114106e-20 * x^9 - 4.267890251195303e-17 * x^8 + 2.879118523083018e-14 * x^7 ...
                - 1.084418008735459e-11 * x^6 + 2.491796224784373e-09 * x^5 - 3.591462823163767e-07 * x^4 ...
                + 3.232810400304542e-05 * x^3 - 1.584729282376364e-03 * x^2 + 5.228413840446568e-02 * x ...
                - 6.547482267336220e-01;
            
            % use formula deducted from Bragg Kleeman rule to calcuate
            % energy straggling given the total sigma and the range
            % straggling
            energySpread = (totalSigmaSq - totalSpreadSq(r80)) / (0.022^2 * 1.77^2 * mcDataEnergy.MeanEnergy^(2*1.77-2));
            energySpread(energySpread < 0) = 0;
            mcDataEnergy.EnergySpread = sqrt(energySpread);
        end
        
        function mcDataOptics = fitBeamOpticsForEnergy(obj,energyIx, focusIndex)
            %function to calculate beam optics used by mcSquare for given
            %energy
            
            i = energyIx;
            
            %calculate geometric distances and extrapolate spot size at nozzle
            SAD = obj.machine.meta.SAD;
            z     = -(obj.machine.data(i).initFocus.dist(focusIndex,:) - SAD);
            sigmaSq = obj.machine.data(i).initFocus.sigma(focusIndex,:).^2;
            
            %fit Courant-Synder equation to data using ipopt, formulae
            %given in mcSquare documentation
            sigmaNull = sqrt(interp1(z,sigmaSq,0));
            
            qRes = @(rho, sigmaT) (sigmaSq -  (sigmaNull^2 - 2*sigmaNull*rho*sigmaT.*z + sigmaT^2.*z.^2));
            
            funcs.objective = @(x) sum(qRes(x(1), x(2)).^2);
            funcs.gradient  = @(x) [  2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(2) * z));
                2 * sum(qRes(x(1), x(2)) .* (2 * sigmaNull * x(1) * z  - 2 * x(2) * z.^2))];
            
            options.lb = [-0.99, -Inf];
            options.ub = [ 0.99,  Inf];
            
            options.ipopt.hessian_approximation = 'limited-memory';
            options.ipopt.limited_memory_update_type = 'bfgs';
            options.ipopt.print_level = 1;
            
            start = [0.9; 0.1];
            [result, ~] = ipopt (start, funcs, options);
            rho    = result(1);
            sigmaT = result(2);
            
            %calculate divergence, spotsize and correlation at nozzle
            DivergenceAtNozzle  = sigmaT;
            SpotsizeAtNozzle    = sqrt(sigmaNull^2 - 2 * rho * sigmaNull * sigmaT * obj.nozzleToIso + sigmaT^2 * obj.nozzleToIso^2);
            CorrelationAtNozzle = (rho * sigmaNull - sigmaT * obj.nozzleToIso) / SpotsizeAtNozzle;
            
            
            %save calcuated beam optics data in mcData
            mcDataOptics.ProtonsMU     = 1e6;
            
            mcDataOptics.Weight1       = 1;
            mcDataOptics.SpotSize1x    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1x  = DivergenceAtNozzle;
            mcDataOptics.Correlation1x = CorrelationAtNozzle;
            mcDataOptics.SpotSize1y    = SpotsizeAtNozzle;
            mcDataOptics.Divergence1y  = DivergenceAtNozzle;
            mcDataOptics.Correlation1y = CorrelationAtNozzle;
            
            mcDataOptics.Weight2       = 0;
            mcDataOptics.SpotSize2x    = 0;
            mcDataOptics.Divergence2x  = 0;
            mcDataOptics.Correlation2x = 0;
            mcDataOptics.SpotSize2y    = 0;
            mcDataOptics.Divergence2y  = 0;
            mcDataOptics.Correlation2y = 0;
            mcDataOptics.FWHMatIso = 2.355 * sigmaNull;
        end
        
        function obj = writeTopasData(obj,ct,stf,pln,w,TopasConfig)
            %function that writes a data file containing stf specific data
            %for a Monte Carlo simulation with TOPAS
            
            %set default values for the simulation, if not specified
            if nargin < 6
                TopasConfig = struct;
            end
            if ~isfield(TopasConfig,'filepath')
                % default: 'topas/MCexport/'
                TopasConfig.filepath = 'topas/MCexport/';
            end
            
            if ~isfield(TopasConfig,'fracHistories')
                % default: 1e-4
                TopasConfig.fracHistories = 1e-4;
            end
            
            if ~isfield(TopasConfig,'minRelWeight')
                % default: 0.001
                % minRelWeight = 0 means all weights are being considered
                % can otherwise be assigned to min(w)
                TopasConfig.minRelWeight = 0.001;
            end
            
            if ~isfield(TopasConfig,'numOfRuns')
                % default: 5 runs
                TopasConfig.numOfRuns = 5;
            end
            
            if ~isfield(TopasConfig,'useOrigBaseData')
                % base data of the original matRad plan will be used
                % default: true
                TopasConfig.useOrigBaseData = false;
            end
            
            if ~isfield(TopasConfig,'beamProfile')
                % default: 5 runs
                TopasConfig.beamProfile = 'biGaussian';
            end
            
            if ~isfield(TopasConfig,'pencilBeamScanning')
                %   beamProfile: 'simple' of 'biGaussian'
                % default: true
                TopasConfig.pencilBeamScanning = 'true';
            end
            
            if ~isfield(TopasConfig,'electronProdCut')
                % default: 0.5
                TopasConfig.electronProdCut = '0.5';
            end
            
            % Check if contradicting settings are being used
            if TopasConfig.useOrigBaseData == true && ~strcmp(TopasConfig.beamProfile,'simple')
                error('Original base data only usable with simple beam geometry!')
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
                
                if TopasConfig.useOrigBaseData
                    [~,ixTmp,~] = intersect([ obj.machine.data.energy], [stf.ray.energy]);
                    for i = 1:length(ixTmp)
                        selectedData(i) =  obj.machine.data(ixTmp(i));
                    end
                    energies = [selectedData.energy];
                else
                    selectedData = [];
                    for i = 1:numel(focusIndex)
                        selectedData = [selectedData, obj.mcSquareData(focusIndex(i), i)];
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
                        minNbParticlesSpot = round(max([TopasConfig.minRelWeight*maxParticlesInSpot,1]));
                        
                        % check whether there are (enough) particles for beam delivery
                        if (voxel_nbParticles>minNbParticlesSpot)
                            
                            cutNumOfBixel = cutNumOfBixel + 1;
                            bixelEnergy = stf(beamIx).ray(rayIx).energy(bixelIx);
                            [~,ixTmp,~] = intersect(energies, bixelEnergy);
                            
                            
                            voxel_x = -stf(beamIx).ray(rayIx).rayPos_bev(3);
                            voxel_y = stf(beamIx).ray(rayIx).rayPos_bev(1);
                            
                            dataTOPAS(cutNumOfBixel).posX = -1.*voxel_x;
                            dataTOPAS(cutNumOfBixel).posY = voxel_y;
                            
                            dataTOPAS(cutNumOfBixel).current = uint32(TopasConfig.fracHistories*voxel_nbParticles);
                            
                            if TopasConfig.pencilBeamScanning
                                % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
                                % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
                                % note that Y' corresponds to the Y axis after the rotation of angleX around X axis
                                dataTOPAS(cutNumOfBixel).angleX = atan(dataTOPAS(cutNumOfBixel).posY / SAD_mm);
                                dataTOPAS(cutNumOfBixel).angleY = atan(-dataTOPAS(cutNumOfBixel).posX ./ (SAD_mm ./ cos(dataTOPAS(cutNumOfBixel).angleX)));
                                dataTOPAS(cutNumOfBixel).posX = (dataTOPAS(cutNumOfBixel).posX / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);
                                dataTOPAS(cutNumOfBixel).posY = (dataTOPAS(cutNumOfBixel).posY / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);
                            end
                            
                            if TopasConfig.useOrigBaseData
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
                
                historyCount = uint32(TopasConfig.fracHistories * nbParticlesTotal);

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
                fileID = fopen([TopasConfig.filepath,'beamSetup_matRad_plan_field',num2str(beamIx),'.txt'],'w');
                
                fprintf(fileID,'i:Ts/ShowHistoryCountAtInterval = %i\n',historyCount/20);
                fprintf(fileID,'s:Sim/PlanLabel = "simData_matRad_plan_field1_run" + Ts/Seed\n');
                fprintf(fileID,'d:Sim/GantryAngle = %.6f deg\n', stf(beamIx).gantryAngle);
                fprintf(fileID,'d:Sim/CouchAngle = %.6f deg\n', stf(beamIx).couchAngle);
                fprintf(fileID,'s:Sim/ParticleName = "proton"\n');
                fprintf(fileID,'u:Sim/ParticleMass = 1.0\n');
                fprintf(fileID,'i:Sim/NbThreads = 0\n');
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
                
                switch TopasConfig.beamProfile
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
                
                if TopasConfig.pencilBeamScanning
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
                if TopasConfig.pencilBeamScanning
                    fprintf(fileID,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
                    fprintf(fileID,'d:Ge/Nozzle/RotZ = 0.0 rad\n');
                end
                
                fprintf(fileID,['d:Ph/Default/CutForElectron = ',TopasConfig.electronProdCut,' mm\n']);
                
                switch TopasConfig.beamProfile
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
                fclose(fileID);
                
                %write run scripts for TOPAS
                basematerial = '';
                if ~exist('machine') || ~isfield( obj.machine.meta,'basematerial')
                    disp('Using default base material: Water')
                    basematerial = 'Water';
                else
                    basematerial =  obj.machine.meta.basematerial;
                end
                
                rspCube = matRad_exportCtTOPAS(ct, TopasConfig.filepath, basematerial);
                
                for runIx = 1:TopasConfig.numOfRuns
                    fileID = fopen([TopasConfig.filepath,'matRad_plan_field',num2str(beamIx),'_run',num2str(runIx),'.txt'],'w');
                    fprintf(fileID,['i:Ts/Seed = ',num2str(runIx),'\n']);
                    fprintf(fileID,'includeFile = ./beamSetup_matRad_plan_field1.txt');
                    fclose(fileID);
                end
                
                %write MCparam file with basic parameters
                MCparam.nbRuns = TopasConfig.numOfRuns;
                MCparam.nbFields = length(stf);
                MCparam.tallies = {'physicalDose'};
                MCparam.simLabel = 'matRad_plan';
                MCparam.cubeDim = ct.cubeDim;
                MCparam.RSP = rspCube;
                MCparam.voxelDimensions = ct.resolution;
                MCparam.nbHistories = historyCount;
                MCparam.nbParticles = nbParticlesTotal;
                save([TopasConfig.filepath,'MCparam.mat'],'MCparam');
            end
        end
        
        function obj = writeMCsquareData(obj,filepath)
            %function that writes a data file containing Monte Carlo base
            %data for a simulation with MCsquare
            
            %look up focus indices
            focusIndex = obj.selectedFocus(obj.energyIndex);
            
            %save mcData acording to used focus index in selectedData
            selectedData = [];
            for i = 1:numel(focusIndex)
                
                selectedData = [selectedData, obj.mcSquareData(focusIndex(i), i)];
            end
            
            %remove field not needed for MCsquare base data
            selectedData = rmfield(selectedData, 'FWHMatIso');
            
            %write MCsqaure data base file
            try
                
                fileID = fopen(filepath,'w');
                
                %Header
                %fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
                fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
                fprintf(fileID,'# %s\n', obj.machine.meta.description);
                fprintf(fileID,'# created by %s on %s\n\n', obj.machine.meta.created_by, obj.machine.meta.created_on);
                
                fprintf(fileID,'Nozzle exit to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.nozzleToIso);
                
                fprintf(fileID,'SMX to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smx);
                
                fprintf(fileID,'SMY to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n',obj.smy);
                
                fprintf(fileID,'Beam parameters\n%d energies\n\n',size(selectedData,2));
                
                fn = fieldnames(selectedData);
                for names = 1:size(fn,1)
                    fprintf(fileID, fn{names});
                    fprintf(fileID, '\t');
                end
                fprintf(fileID, '\n');
                
                for k = 1:size(selectedData,2)
                    for m = 1:numel(fn)
                        fprintf(fileID, '%g', selectedData(k).(fn{m}));
                        fprintf(fileID, '\t');
                    end
                    fprintf(fileID, '\n');
                end
                
                fclose(fileID);
                
                obj.bdl_path = filepath;
                
            catch MException
                error(MException.message);
            end
        end
        
        function obj = saveMatradMachine(obj,name)
            %save previously calculated mcSquareData in new baseData file
            %with given name
            
            [~ ,energyIndex, ~] = intersect([obj.machine.data(:).energy], [obj.mcSquareData(:).NominalEnergy]);
            
            machineName = [obj.machine.meta.radiationMode, '_', name];
            
            count = 1;
            for i = energyIndex'
                
                obj.machine.data(i).mcSquareData = obj.mcSquareData(:,count);
                
                count = count + 1;
            end
            
            save(strcat('../../', machineName, '.mat'),'machine');
        end
    end
    
end

