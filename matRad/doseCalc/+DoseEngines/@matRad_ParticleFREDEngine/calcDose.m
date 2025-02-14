function dij = calcDose(this,ct,cst,stf)
% Function to forward dose calculation to FRED and inport the results
% in matRad
%
% call
%   dij = this.calcDose(ct,stf,pln,cst)
%
% input
%   ct:          	matRad ct struct
%   cst:            matRad cst struct
%   stf:         	atRad steering information struct
%   
% output
%   dij:            matRad dij struct
%
% References
%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    matRad_cfg = MatRad_Config.instance();
    
    currFolder = pwd;
    cd(this.FREDrootFolder);

    %Now we can run initDoseCalc as usual
    dij = this.initDoseCalc(ct,cst,stf);

    % Interpolate cube on dose grid
    HUcube{1} =  matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y',  dij.ctGrid.z,ct.cubeHU{1}, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'linear');

    % Force HU clamping if values are found outside of available range
    switch this.HUtable
       case 'internal'
           if any(HUcube{1}(:)>this.hLutLimits(2)) || any(HUcube{1}(:)<this.hLutLimits(1))
               matRad_cfg.dispWarning('HU outside of boundaries');
               this.HUclamping = true;
           end
       otherwise
           matRad_cfg.dispInfo('Using custom HU table: %s\n', this.HUtable);
    end
    
    if this.ignoreOutsideDensities
        eraseCtDensMask = ones(prod(ct.cubeDim),1);
        eraseCtDensMask(this.VctGrid) = 0;
        HUcube{1}(eraseCtDensMask == 1) = this.hLutLimits(1);
    end

    %Write the directory tree
    this.writeTreeDirectory();
    
    % Write patient
    cd(this.regionsFolder);

    fileNamePatient = fullfile(this.regionsFolder, this.patientFilename);
    patientMetadata.imageOrigin = [0 0 0];
    patientMetadata.resolution  = [this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z];

    patientMetadata.datatype = 'int16';

    %%%%%% !!!!!!!!!!!!! mind this flip !!!!!!!!!!!!! %%%%%
    % Need to permute x and y because of the order data is written in mhd
    % HUcube{1} = permute(HUcube{1}, [2,1,3]);
    matRad_writeMHD(fileNamePatient, HUcube{1},patientMetadata);

    cd(this.FREDrootFolder);

    % Linear projection of BEV source (x,y) points to plane at BAMStoISO distance 
    getPointAtBAMS = @(target,source,distance,BAMStoIso) (target -source)*(-BAMStoIso)/distance + source;

    % Loop over the stf to rearrange data
    counter = 0;

    % Use the emittance base data class to recover MC information
    emittanceBaseData = matRad_MCemittanceBaseData(this.machine,stf);
    
    
    % Loop over fields. FRED performs one single simulation for multiple
    % fields
    for i = 1:length(stf)

        stfFred(i).gantryAngle     = stf(i).gantryAngle;
        stfFred(i).couchAngle      = stf(i).couchAngle;
        
        % Get dose grid resolution
        doseGridResolution = [this.doseGrid.resolution.x, this.doseGrid.resolution.y, this.doseGrid.resolution.z];

        % get matRad isocenter coordinates in dose cube coordinate system.
        % -> Dose cube coordinate system places the first coordinate point (aka resolution) in
        % the center of the first voxel. So the "zero" of the cube
        % coordinate system is 0.5*resolution outside of the cube surface.
        isoInDoseGridCoord = matRad_world2cubeCoords(stf(i).isoCenter,this.doseGrid);

        % Coordinate of the first voxel in cube system is in the center of the first voxel and is equal to resolution. Thus the zero 
        % of the cube coordinate system is 0.5*resolution before the surface of
        % the phantom. The surface thus starts at 1/2 resolution wrt zero of that system
        fredCubeSurfaceInDoseCubeCoords =  0.5*doseGridResolution;
        
        % Get coordinates of pivot point in FRED cube (center of
        % geometrical cube) in dose cube coordinates. This is the distance
        % between pivot point and the surface + the position of the
        % surface in the cube coord. system.
        fredPivotInCubeCoordinates = 0.5*this.doseGrid.dimensions([2 1 3]).*doseGridResolution + fredCubeSurfaceInDoseCubeCoords;
        
        
        % Define the FRED isocenter as the distance between the pivot point
        % and the matRad isocenter in the cube coordinate system.
        stfFred(i).isoCenter = -(fredPivotInCubeCoordinates - isoInDoseGridCoord);

        % First coordinate is flipped
        stfFred(i).isoCenter = stfFred(i).isoCenter.*[-1 1 1];
        
        % NOTE on the coordinate system.
        % FRED places the pivot point of the component at the center of the
        % FRED coordinate system, then applies a translation s.t. the FRED
        % isocenter (defined for each field) is in the center of the FRED
        % coordinate system. Then applies rotations. This way everything
        % is defined in BEV reference.

        nominalEnergies        = unique([stf(i).ray.energy]);
        [~,nominalEnergiesIdx] = intersect([this.machine.data.energy],nominalEnergies);
        
        energyIdxInEmittance   = ismember(emittanceBaseData.energyIndex, nominalEnergiesIdx);
        monteCarloBaseData     = emittanceBaseData.monteCarloData(energyIdxInEmittance);
        
        stfFred(i).nominalEnergies = nominalEnergies;
        stfFred(i).energies        = [monteCarloBaseData.MeanEnergy].*this.numOfNucleons.*this.primaryMass;        % Note for generic baseData: the kernels were simulated with equivalent of primaryMass = 1
        stfFred(i).energySpread    = [monteCarloBaseData.EnergySpread];
        stfFred(i).energySpreadMeV = [monteCarloBaseData.EnergySpread].*[monteCarloBaseData.MeanEnergy]/100;
        stfFred(i).FWHMs = 2.355*[monteCarloBaseData.SpotSize1x];
        
        stfFred(i).energySpreadFWHMMev    = 2.355*stfFred(i).energySpreadMeV;
        stfFred(i).BAMStoIsoDist   = emittanceBaseData.nozzleToIso;
        
        % Select the parametrs for source model
        switch this.sourceModel

            case 'gaussian'
        
            case 'emittance'
                stfFred(i).emittanceX = [];
                stfFred(i).twissBetaX  = [];
                stfFred(i).twissAlphaX = [];
                stfFred(i).emittanceRefPlaneDistance = [];
                
                % Need to get the parameters for the model from MCemittance
                for eIdx=emittanceBaseData.energyIndex'
                    % Only using first focus index for now
                    tmpOpticsData = emittanceBaseData.fitBeamOpticsForEnergy(eIdx,1);
                    stfFred(i).emittanceX       = [stfFred(i).emittanceX,  tmpOpticsData.twissEpsilonX];
                    stfFred(i).twissBetaX       = [stfFred(i).twissBetaX,  tmpOpticsData.twissBetaX];
                    stfFred(i).twissAlphaX      = [stfFred(i).twissAlphaX, tmpOpticsData.twissAlphaX];
                    stfFred(i).emittanceRefPlaneDistance = [stfFred(i).emittanceRefPlaneDistance, this.machine.meta.BAMStoIsoDist];
                end

            case 'sigmaSqrModel'
                stfFred(i).sSQr_a = [];
                stfFred(i).sSQr_b = [];
                stfFred(i).sSQr_c = [];

                for eIdx=emittanceBaseData.energyIndex'
                    tmpOpticsData = emittanceBaseData.fitBeamOpticsForEnergy(eIdx,1);
                    stfFred(i).sSQr_a           = [stfFred(i).sSQr_a, tmpOpticsData.sSQ_a];
                    stfFred(i).sSQr_b           = [stfFred(i).sSQr_b, tmpOpticsData.sSQ_b];
                    stfFred(i).sSQr_c           = [stfFred(i).sSQr_c, tmpOpticsData.sSQ_c];
                end
            otherwise
                matRad_cfg.dispWarning('Unrecognized source model, setting gaussian');

        end

        % Allocate empty layer container
        % Rearrange info into separate energy layers
        for j = 1:numel(stfFred(i).energies)

            %stfFred(i).energyLayer(j).targetPoints   = [];
            stfFred(i).energyLayer(j).numOfPrimaries = [];
            stfFred(i).energyLayer(j).rayNum         = [];
            stfFred(i).energyLayer(j).bixelNum       = [];
            stfFred(i).energyLayer(j).rayDivX        = [];
            stfFred(i).energyLayer(j).rayDivY        = [];
            stfFred(i).energyLayer(j).rayPosX        = [];
            stfFred(i).energyLayer(j).rayPosY        = [];
        end

        for j = 1:stf(i).numOfRays
            for k = 1:stf(i).numOfBixelsPerRay(j)
                counter = counter + 1;
                dij.beamNum(counter,1)  = i;
                dij.rayNum(counter,1)   = j;
                dij.bixelNum(counter,1) = k;
            end
           
            for k = 1:numel(stfFred(i).energies)

                if any(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))
                    stfFred(i).energyLayer(k).rayNum   = [stfFred(i).energyLayer(k).rayNum j];

                    stfFred(i).energyLayer(k).bixelNum = [stfFred(i).energyLayer(k).bixelNum ...
                        find(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))];
                    
                    % Get spot position and divergence
                    targetX = stf(i).ray(j).targetPoint_bev(1);
                    targetY = stf(i).ray(j).targetPoint_bev(3);

                    % Stf.ray.rayPos_bev is position of ray at the
                    % IsoCenter plane.
                    sourceX = stf(i).ray(j).rayPos_bev(1);
                    sourceY = stf(i).ray(j).rayPos_bev(3);

                    distance = stf(i).ray(j).targetPoint_bev(2) - stf(i).ray(j).rayPos_bev(2);
                    
                    divergenceX = (targetX - sourceX)/distance;
                    divergenceY = (targetY - sourceY)/distance;
                    
                    %stfFred(i).energyLayer(k).targetPoints    = [stfFred(i).energyLayer(k).targetPoints; -targetX targetY];

                    % This is position of the spot at -BAMsToIso distance
                    % (zero is at IsoCenter depth).
                    stfFred(i).energyLayer(k).rayPosX         = [stfFred(i).energyLayer(k).rayPosX, getPointAtBAMS(targetX,sourceX,distance,stfFred(i).BAMStoIsoDist)];
                    stfFred(i).energyLayer(k).rayPosY         = [stfFred(i).energyLayer(k).rayPosY, getPointAtBAMS(targetY,sourceY,distance,stfFred(i).BAMStoIsoDist)];

                    stfFred(i).energyLayer(k).rayDivX         = [stfFred(i).energyLayer(k).rayDivX, divergenceX];
                    stfFred(i).energyLayer(k).rayDivY         = [stfFred(i).energyLayer(k).rayDivY, divergenceY];
                    
                    
                    if this.calcDoseDirect
                        % Set the bixel weight
                        stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries ...
                                                  stf(i).ray(j).weight(stf(i).ray(j).energy == stfFred(i).nominalEnergies(k))];
                    else
                         stfFred(i).energyLayer(k).numOfPrimaries = [stfFred(i).energyLayer(k).numOfPrimaries, 1];
                    end

                end

            end
        end

        %FRED works in cm
        stfFred(i).isoCenter       = stfFred(i).isoCenter/10;
        stfFred(i).BAMStoIsoDist   = stfFred(i).BAMStoIsoDist/10;

        switch this.sourceModel
            case 'gaussian'
                stfFred(i).FWHMs           = stfFred(i).FWHMs/10;                
            case 'emittance'
                stfFred(i).emittanceRefPlaneDistance = stfFred(i).emittanceRefPlaneDistance/10;
            case 'sigmaSqrModel'

        end

        stfFred(i).totalNumOfBixels = stf(i).totalNumOfBixels;
        for j=1:numel(stfFred(i).nominalEnergies)
           stfFred(i).energyLayer(j).rayPosX      = stfFred(i).energyLayer(j).rayPosX/10;
           stfFred(i).energyLayer(j).rayPosY      = stfFred(i).energyLayer(j).rayPosY/10;
           stfFred(i).energyLayer(j).nBixels      = numel(stfFred(i).energyLayer(j).bixelNum);

           if this.calcDoseDirect
               stfFred(i).energyLayer(j).numOfPrimaries = this.conversionFactor*stfFred(i).energyLayer(j).numOfPrimaries;
           end
        end
    end
    
    counterFred = 0;
    fredOrder = NaN * ones(dij.totalNumOfBixels,1);
    for i = 1:length(stf)
        for j = 1:numel(stfFred(i).nominalEnergies)
            for k = 1:numel(stfFred(i).energyLayer(j).numOfPrimaries)
                counterFred = counterFred + 1;
                ix = find(i                               == dij.beamNum & ...
                    stfFred(i).energyLayer(j).rayNum(k)   == dij.rayNum & ...
                    stfFred(i).energyLayer(j).bixelNum(k) == dij.bixelNum);
    
                fredOrder(ix) = counterFred;
            end
        end
    end
    
    if any(isnan(fredOrder))
        matRad_cfg.dispError('Invalid ordering of Beamlets for FRED computation!');
    end
    
    % %% MC computation and dij filling
    this.writeFredInputAllFiles(stfFred);

    switch this.externalCalculation
        
        case 'write' % Write simulation files for external calculation (no FRED installation required)
    
            matRad_cfg.dispInfo('All files have been generated\n');
            dijFieldsToOverride = {'numOfBeams','beamNum','bixelNum','rayNum','totalNumOfBixels','totalNumOfRays','numOfRaysPerBeam'};
                
            for fieldName=dijFieldsToOverride
                dij.(fieldName{1}) = this.numOfColumnsDij;
            end

            doseCube = [];
        
        case 'off' % Run FRED simulation (requires installation)

            % Check consistency of installation
            if this.checkExec()
                
                matRad_cfg.dispInfo('calling FRED');
    
                cd(this.MCrunFolder);
                
                systemCall = [this.cmdCall, '-f fred.inp'];
                if ~this.useGPU
                    systemCall = [this.cmdCall, ' -nogpu -f fred.inp'];
                end
    
                % printOutput to matlab console
                if this.printOutput
                    [status,~] = system(systemCall,'-echo');
                else
                    [status,~] = system(systemCall);
                end
                cd(this.FREDrootFolder);
            else
                matRad_cfg.dispError('FRED setup incorrect for this plan simulation');
            end

            if status==0
                matRad_cfg.dispInfo(' done\n');
            end
            
            % read simulation output
            [doseCube, letdCube] = this.readSimulationOutput(this.MCrunFolder,this.calcDoseDirect, 'calcLET', logical(this.calcLET), 'readFunctionHandle', this.dijReaderHandle);

        otherwise % A path for loading has been provided
            
            matRad_cfg.dispInfo(['Reading simulation data from: ', strrep(this.MCrunFolder,'\','\\'), '\n']);

            % read simulation output
            [doseCube, letdCube, loadFileName] = this.readSimulationOutput(this.MCrunFolder,this.calcDoseDirect, 'calcLET',logical(this.calcLET),'readFunctionHandle', this.dijReaderHandle);

            dij.externalCalculationLodPath = loadFileName;

    end

    if ~isempty(doseCube)
    
        % Fill dij
        if this.calcDoseDirect
            % Dose cube
            if isequal(size(doseCube), this.doseGrid.dimensions)
                dij.physicalDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),doseCube(this.VdoseGrid), this.doseGrid.numOfVoxels,1);
            end
    
            % LETd cube
            if this.calcLET
                if isequal(size(letdCube), this.doseGrid.dimensions)
                    dij.mLETd{1}    = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),letdCube(this.VdoseGrid)./10, this.doseGrid.numOfVoxels,1);

                    % We need LETd * dose as well
                    dij.mLETDose{1} = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),(letdCube(this.VdoseGrid)./10).*doseCube(this.VdoseGrid), this.doseGrid.numOfVoxels,1);
                end
            end
    
            % Needed for calcCubes
            dijFieldsToOverride = {'numOfBeams','beamNum','bixelNum','rayNum','totalNumOfBixels','totalNumOfRays','numOfRaysPerBeam'};
            
            for fieldName=dijFieldsToOverride
                dij.(fieldName{1}) = 1;
            end

        else
            % Dose cube
            if isequal(size(doseCube), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                %When scoring dij, FRED internaly normalizes to 1
                dij.physicalDose{1}(this.VdoseGrid,:) = this.conversionFactor*doseCube(this.VdoseGrid,fredOrder);
            end

            % LET cube
            if this.calcLET
                if isequal(size(letdCube), [dij.doseGrid.numOfVoxels,dij.totalNumOfBixels])
                    % Need to divide by 10, FRED scores in MeV * cm^2 / g
                    dij.mLETd{1}(this.VdoseGrid,:) = letdCube(this.VdoseGrid,fredOrder)./10;
                end

                % We need LETd * dose as well
                dij.mLETDose{1} = sparse(dij.physicalDose{1}.*dij.mLETd{1});
            end
        end


        % Calc Biological quantities
        if this.calcBioDose
           % recover alpha and beta maps
           tmpBixel.radDepths = zeros(size(this.VdoseGrid,1),1);
        
           tmpBixel.vAlphaX   = dij.ax{1}(this.VdoseGrid);
           tmpBixel.vBetaX    = dij.bx{1}(this.VdoseGrid);
           tmpBixel.vABratio  = dij.ax{1}(this.VdoseGrid)./dij.bx{1}(this.VdoseGrid);

           if this.calcDoseDirect
                tmpKernel.LET = dij.mLETd{1}(this.VdoseGrid);

                tmpBixel = this.bioModel.calcBiologicalQuantitiesForBixel(tmpBixel,tmpKernel);
                
                tmpBixel.alpha(isnan(tmpBixel.alpha)) = 0;
                tmpBixel.beta(isnan(tmpBixel.beta)) =  0;

                dij.mAlphaDose{1}     = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),tmpBixel.alpha.*dij.physicalDose{1}(this.VdoseGrid), this.doseGrid.numOfVoxels,1);
                dij.mSqrtBetaDose{1}  = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),sqrt(tmpBixel.beta).*dij.physicalDose{1}(this.VdoseGrid), this.doseGrid.numOfVoxels,1);
           else    
               % Loop over all bixels
               for bxlIdx = 1:dij.totalNumOfBixels
                   bixelLET           = full(dij.mLETd{1}(:,bxlIdx));
                   tmpKernel.LET      = bixelLET(this.VdoseGrid);

                   tmpBixel = this.bioModel.calcBiologicalQuantitiesForBixel(tmpBixel,tmpKernel);

                   tmpBixel.alpha(isnan(tmpBixel.alpha)) = 0;
                   tmpBixel.beta(isnan(tmpBixel.beta)) =  0;
                
                   dij.mAlphaDose{1}(:,bxlIdx)     = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),tmpBixel.alpha.*dij.physicalDose{1}(this.VdoseGrid,bxlIdx), this.doseGrid.numOfVoxels,1);
                   dij.mSqrtBetaDose{1}(:,bxlIdx)  = sparse(this.VdoseGrid, ones(numel(this.VdoseGrid),1),sqrt(tmpBixel.beta).*dij.physicalDose{1}(this.VdoseGrid,bxlIdx), this.doseGrid.numOfVoxels,1);
               end
           end
        end
    end

    dij = this.finalizeDose(dij);
    
    cd(currFolder);
end