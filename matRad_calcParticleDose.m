function dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
% matRad particle dose calculation wrapper
%
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst,calcDoseDirect)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: (optional) switches to direct dose calculation, only
%                   makes sense in combination with matRad_calcDoseDirect
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% initialize waitbar
if matRad_cfg.logLevel > 1
    figureWait = waitbar(0,'calculate dose influence matrix for particles...');
    % prevent closure of waitbar and show busy state
    set(figureWait,'pointer','watch');
end

matRad_cfg.dispInfo('matRad: Particle dose calculation... \n');

% load default parameters in case they haven't been set yet
pln = matRad_cfg.getDefaultProperties(pln,{'propDoseCalc'});

if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
    pln.propHeterogeneity = matRad_HeterogeneityConfig();
    pln.propHeterogeneity.bioOpt = pln.bioParam.bioOpt;
    matRad_cfg.dispInfo(['Modulation power set to Pmod = ' num2str(pln.propHeterogeneity.modPower) ' µm.\n']);
    cstOriginal = cst;
end

% init dose calc
matRad_calcDoseInit;

% initialize lung heterogeneity correction and turn off if necessary files are missing
if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
    matRad_cfg.dispInfo('Heterogeneity correction enabled. \n');
    heteroCST = false;
    for i = 1:length(cst(:,1)) % scan cst for segmentation flagged for correction
        if isfield(cst{i,5},'HeterogeneityCorrection')
            heteroCST = true;
            continue
        end
    end
    if ~isstruct(machine.data(1).Z) || ~heteroCST
        matRad_cfg.dispWarning('Heterogeneity correction enabled but no usable data in cst or unsuitable base data. Correction cannot be applied.');
        pln.propHeterogeneity.calcHetero = false;
    end
else
    matRad_cfg.dispInfo('Heterogeneity correction disabled. \n');
end

% initialize HeteroCorrStruct and adjust base data if needed
if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
    if pln.bioParam.bioOpt && ~isfield(machine.data,'alpha')
        matRad_cfg.dispInfo('Calculating alpha-beta curves for baseData ... ');
        machine = matRad_getAlphaBetaCurves(machine,pln,cst);
        matRad_cfg.dispInfo('Done!\n');
    end

    % get all lung voxel indices
    lungVoxel = unique(cell2mat([cstOriginal{cellfun(@(teststr) ~isempty(strfind(lower(teststr),'lung')), cst(:,2)),4}]'),'rows');
    calcHeteroCorrStruct.cubeDim = ct.cubeDim;
    calcHeteroCorrStruct.numOfCtScen = pln.multScen.numOfCtScen;
    calcHeteroCorrStruct.resolution = ct.resolution;

    calcHeteroCorrStruct.cube = cell(1,ctScen);
    calcHeteroCorrStruct.cube(1,:) = {zeros(ct.cubeDim)};

    for shiftScen = 1:pln.multScen.numOfCtScen
        calcHeteroCorrStruct.cube{shiftScen}(lungVoxel(:,shiftScen)) = ct.cube{shiftScen}(lungVoxel(:,shiftScen));
    end
end

if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero && pln.propHeterogeneity.useOriginalDepths
    machine.data = matRad_HeterogeneityConfig.overrideBaseData(machine.data);
end

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% if biological optimization considering a variable RBE is true then create alphaDose and betaDose containers and sparse matrices
if pln.bioParam.bioOpt

    alphaDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);

    for ctScen = 1:pln.multScen.numOfCtScen
        for shiftScen = 1:pln.multScen.totNumShiftScen
            for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                    dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}        = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                    dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}     = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                end
            end
        end
    end
end

if pln.propDoseCalc.calcLET
    if isfield(machine.data,'LET')

        letDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);

        for ctScen = 1:pln.multScen.numOfCtScen
            for shiftScen = 1:pln.multScen.totNumShiftScen
                for rangeShiftScen = 1:pln.multScen.totNumRangeScen

                    if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                        dij.mLETDose{ctScen,shiftScen,rangeShiftScen} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
                    end

                end
            end
        end
        matRad_cfg.dispInfo('LET computation enabled!\n');
    else
        matRad_cfg.dispWarning('LET not available in the machine data. LET will not be calculated.');
    end
end

%Toggles correction of small difference of current SSD to distance used
%in generation of base data (e.g. phantom surface at isocenter)
if pln.propDoseCalc.airOffsetCorrection
    if ~isfield(machine.meta, 'fitAirOffset')
        fitAirOffset = 0; %By default we assume that the base data was fitted to a phantom with surface at isocenter
        matRad_cfg.dispDebug('Asked for correction of Base Data Air Offset, but no value found. Using default value of %f mm.\n',fitAirOffset);
    else
        fitAirOffset = machine.meta.fitAirOffset;
    end
else
    fitAirOffset = 0;
end

if ~isfield(machine.meta, 'BAMStoIsoDist')
    BAMStoIsoDist = 1000;
    matRad_cfg.dispWarning('Machine data does not contain BAMStoIsoDist. Using default value of %f mm.',BAMStoIsoDist);
else
    BAMStoIsoDist = machine.meta.BAMStoIsoDist;
end

% book keeping - this is necessary since pln is not used in optimization or
% matRad_calcCubes
if strcmp(pln.bioParam.model,'constRBE')
    dij.RBE = pln.bioParam.RBE;
end


% generates tissue class matrix for biological treatment planning and alpha_x, beta_x, vectors
if pln.bioParam.bioOpt

    vTissueIndex = zeros(size(VdoseGrid,1),1);
    dij.ax       = zeros(dij.doseGrid.numOfVoxels,ct.numOfCtScen);
    dij.bx       = zeros(dij.doseGrid.numOfVoxels,ct.numOfCtScen);
    dij.abx      = zeros(dij.doseGrid.numOfVoxels,ct.numOfCtScen);  % alpha beta ratio

    % retrieve photon LQM parameter for the current dose grid voxels
    %     [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,ct.numOfCtScen,VdoseGrid);
    [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,ct.numOfCtScen);


    dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);

    % only if LEM is used corresponding bio data must be available in the base data set
    if strcmp(pln.bioParam.model,'LEM') || (isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero)
        if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')

            matRad_cfg.dispInfo('loading biological base data...');

            for i = 1:size(cst,1)

                % check if cst is compatiable
                if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX')

                    % check if base data contains alphaX and betaX
                    IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                        ismember(machine.data(1).betaX, cst{i,5}.betaX));

                    % check consistency of biological baseData and cst settings
                    if ~isempty(IdxTissue)
                        isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{1});
                        vTissueIndex(isInVdoseGrid) = IdxTissue;
                    else
                        matRad_cfg.dispError('biological base data and cst inconsistent!');
                    end

                else
                    vTissueIndex(row) = 1;
                    matRad_cfg.dispInfo(' tissue type of %s was set to 1...',cst{i,2});
                end
            end

            matRad_cfg.dispInfo('Done!\n');

        else
            matRad_cfg.dispError('base data is incomplete - alphaX and/or betaX is missing');
        end

    else
        % parametrized biological models are based on the LET
        if ~isfield(machine.data,'LET')
            matRad_cfg.dispError('base data is incomplement - LET is missing');
        end
    end %  end is LEM model
end

% lateral cutoff for raytracing and geo calculations
pln.propDoseCalc.effectiveLateralCutOff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;

if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'lateralCutOff')
    pln.propDoseCalc.lateralCutOff = matRad_cfg.propDoseCalc.defaultLateralCutOff;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over all shift scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for shiftScen = 1:pln.multScen.totNumShiftScen

    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(shiftScen,:);
    end

    if pln.multScen.totNumShiftScen > 1
        matRad_cfg.dispInfo('\tShift scenario %d of %d: \n',shiftScen,pln.multScen.totNumShiftScen);
    end

    counter = 0;

    % compute SSDs only for first scenario
    stf = matRad_computeSSD(stf,ct);

    for i = 1:numel(stf) % loop over all beams

        % init beam
        matRad_calcDoseInitBeam;

        % Calculate radiological depth cube for heterogeneity correction
        if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
            matRad_cfg.dispInfo('matRad: calculate radiological depth cube for heterogeneity correction...');
            heteroCorrDepthV = matRad_rayTracing(stf(i),calcHeteroCorrStruct,VctGrid,rot_coordsV,pln.propDoseCalc.effectiveLateralCutOff);

            % HETERO interpolate hetero depth cube to dose grid resolution
            heteroCorrDepthV = matRad_interpRadDepth...
                (ct,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,heteroCorrDepthV);
            matRad_cfg.dispInfo('Done!\n');
        end

        % Determine lateral cutoff
        matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');
        visBoolLateralCutOff = 0;
        machine = matRad_calcLateralParticleCutOff(machine,pln.propDoseCalc.lateralCutOff,stf(i),visBoolLateralCutOff);
        matRad_cfg.dispInfo('Done!\n');

        for j = 1:stf(i).numOfRays % loop over all rays

            if ~isempty(stf(i).ray(j).energy)

                % find index of maximum used energy (round to keV for numerical reasons
                energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);

                maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);

                % calculate initial sigma for all bixel on current ray
                sigmaIniRay = matRad_calcSigmaIni(machine.data,stf(i).ray(j),stf(i).ray(j).SSD);

                if strcmp(pln.propDoseCalc.fineSampling.calcMode, 'fineSampling')
                    % Ray tracing for beam i and ray j with explicit
                    % lateral distances for fine sampling
                    [ix,~,~,~,latDistsX,latDistsZ] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                        stf(i).sourcePoint_bev, ...
                        stf(i).ray(j).targetPoint_bev, ...
                        machine.meta.SAD, ...
                        find(~isnan(radDepthVdoseGrid{1})), ...
                        maxLateralCutoffDoseCalc);

                    % Given the initial sigmas of the sampling ray, this
                    % function provides the weights for the sub-pencil beams,
                    % their positions and their sigma used for dose calculation
                    for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
                        if (pln.propDoseCalc.fineSampling.sigmaSub < sigmaIniRay(k)) && (pln.propDoseCalc.fineSampling.sigmaSub > 0)
                            [finalWeight(:,k), sigmaSub(:,k), posX(:,k), posZ(:,k), numOfSub(:,k)] = ...
                                matRad_calcWeights(sigmaIniRay(k), pln.propDoseCalc.fineSampling.method, pln.propDoseCalc.fineSampling.N, pln.propDoseCalc.fineSampling.sigmaSub);
                        else
                            if (pln.propDoseCalc.fineSampling.sigmaSub < 0)
                                matRad_cfg.dispError('Chosen fine sampling sigma cannot be negative!');
                            elseif (pln.propDoseCalc.fineSampling.sigmaSub > sigmaIniRay(k))
                                matRad_cfg.dispError('Chosen fine sampling sigma is too high for defined plan!');
                            end
                        end
                    end
                else
                    % Ray tracing for beam i and ray j without explicitly
                    % obtaining lateral distances
                    [ix,currRadialDist_sq,~,~,~,~] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                        stf(i).sourcePoint_bev, ...
                        stf(i).ray(j).targetPoint_bev, ...
                        machine.meta.SAD, ...
                        find(~isnan(radDepthVdoseGrid{1})), ...
                        maxLateralCutoffDoseCalc);
                end

                if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
                    heteroCorrDepths = heteroCorrDepthV{1}(ix);
                end

                % just use tissue classes of voxels found by ray tracer
                if pln.bioParam.bioOpt
                    vTissueIndex_j = vTissueIndex(ix,:);
                else
                    vTissueIndex_j = zeros(size(ix));
                end
                for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                    counter       = counter + 1;
                    bixelsPerBeam = bixelsPerBeam + 1;

                    if matRad_cfg.logLevel > 1
                        % Display progress and update text only 200 times
                        if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                            matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                        end

                        % update waitbar only 100 times if it is not closed
                        if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                            waitbar(counter/dij.totalNumOfBixels,figureWait);
                        end
                    end

                    % remember beam and bixel number
                    if ~calcDoseDirect
                        dij.beamNum(counter)  = i;
                        dij.rayNum(counter)   = j;
                        dij.bixelNum(counter) = k;
                        % extract MU data if present (checks for downwards compatability)
                        minMU = 0;
                        if isfield(stf(i).ray(j),'minMU')
                            minMU = stf(i).ray(j).minMU(k);
                        end
                        maxMU = Inf;
                        if isfield(stf(i).ray(j),'maxMU')
                            maxMU = stf(i).ray(j).maxMU(k);
                        end
                        numParticlesPerMU = 1e6;
                        if isfield(stf(i).ray(j),'numParticlesPerMU')
                            numParticlesPerMU = stf(i).ray(j).numParticlesPerMU(k);
                        end
                        dij.minMU(counter,1) = minMU;
                        dij.maxMU(counter,1) = maxMU;
                        dij.numParticlesPerMU(counter,1) = numParticlesPerMU;
                    end

                    % find energy index in base data
                    energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));

                    % Since matRad's ray cast starts at the skin and base data
                    % is generated at soume source to phantom distance
                    % we can explicitly correct for the nozzle to air WEPL in
                    % the current case.
                    if  pln.propDoseCalc.airOffsetCorrection
                        nozzleToSkin = ((stf(i).ray(j).SSD + BAMStoIsoDist) - machine.meta.SAD);
                        dR = 0.0011 * (nozzleToSkin - fitAirOffset);
                    else
                        dR = 0;
                    end

                    % create offset vector to account for additional offsets modelled in the base data and a potential
                    % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                    % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
                    % interpolations in matRad_calcParticleDoseBixel() and avoid extrapolations.
                    offsetRadDepth = machine.data(energyIx).offset - (stf(i).ray(j).rangeShifter(k).eqThickness + dR);

                    % calculate projected coordinates for fine sampling of each beamlet
                    if strcmp(pln.propDoseCalc.fineSampling.calcMode, 'fineSampling')
                        projCoords = matRad_projectOnComponents(VdoseGrid(ix), size(radDepthsMat{1}), stf(i).sourcePoint_bev,...
                            stf(i).ray(j).targetPoint_bev, stf(i).isoCenter,...
                            [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                            -posX(:,k), -posZ(:,k), rotMat_system_T);
                    end
                    % We do now loop over scenarios that alter voxel
                    % values, e.g. range scenarios or ct phases, as we can
                    % vectorize computations more efficiently than when
                    % making this an outer loop
                    for ctScen = 1:pln.multScen.numOfCtScen
                        if any(any(pln.multScen.scenMask(ctScen,:,:))) %We don't need it if no scenario for this ct scenario is relevant
                            % precomputations for fine-sampling
                            if strcmp(pln.propDoseCalc.fineSampling.calcMode, 'fineSampling')
                                % compute radial distances relative to pencil beam
                                % component
                                currRadialDist_sq = reshape(bsxfun(@plus,latDistsX,posX(:,k)'),[],1,numOfSub(k)).^2 + reshape(bsxfun(@plus,latDistsZ,posZ(:,k)'),[],1,numOfSub(k)).^2;

                                % interpolate radiological depths at projected
                                % coordinates
                                radDepths = interp3(radDepthsMat{ctScen},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                                    projCoords(:,2,:)./dij.doseGrid.resolution.y,projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest');
                            else
                                radDepths = radDepthVdoseGrid{ctScen}(ix);
                            end
                        end
                        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                % manipulate radDepthCube for range scenarios
                                if pln.multScen.relRangeShift(rangeShiftScen) ~= 0 || pln.multScen.absRangeShift(rangeShiftScen) ~= 0
                                    currRadDepths = radDepths * (1+pln.multScen.relRangeShift(rangeShiftScen)) +... % rel range shift
                                        pln.multScen.absRangeShift(rangeShiftScen);                                   % absolute range shift
                                    currRadDepths(currRadDepths < 0) = 0;
                                else
                                    currRadDepths = radDepths;
                                end
                                % find depth depended lateral cut off
                                if pln.propDoseCalc.lateralCutOff >= 1
                                    currIx = currRadDepths <= machine.data(energyIx).depths(end) + offsetRadDepth;
                                elseif pln.propDoseCalc.lateralCutOff < 1 && pln.propDoseCalc.lateralCutOff > 0
                                    % perform rough 2D clipping
                                    currIx = currRadDepths <= machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                        currRadialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                                    % peform fine 2D clipping
                                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                            (machine.data(energyIx).LatCutOff.CutOff.^2)', currRadDepths(currIx)) >= currRadialDist_sq(currIx);
                                    end
                                else
                                    matRad_cfg.dispError('Lateral Cut-Off must be a value between 0 and 1!')
                                end

                                % empty bixels may happen during recalculation of error
                                % scenarios -> skip to next bixel
                                if ~any(currIx)
                                    %Create empty container entries for
                                    %this bixel
                                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(dij.doseGrid.numOfVoxels,1);
                                    if isfield(dij,'mLETDose')
                                        letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(dij.doseGrid.numOfVoxels,1);
                                    end
                                    if pln.bioParam.bioOpt
                                        alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                        betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}  = sparse(VdoseGrid(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                                    end

                                    %skip bixel
                                    continue;
                                end

                                % adjust radDepth according to range shifter
                                if  pln.propDoseCalc.airOffsetCorrection
                                    currRadDepths(currIx) = currRadDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness + dR;

                                    %sanity check due to negative corrections
                                    currRadDepths(currRadDepths < 0) = 0;
                                else
                                    currRadDepths(currIx) = currRadDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;
                                end

                                if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
                                    currHeteroCorrDepths = heteroCorrDepths(currIx);
                                end

                                % select correct initial focus sigma squared
                                sigmaIni_sq = sigmaIniRay(k)^2;

                                % consider range shifter for protons if applicable
                                if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                                    % compute!
                                    sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx), ...
                                        stf(i).ray(j).rangeShifter(k), ...
                                        stf(i).ray(j).SSD);

                                    % add to initial sigma in quadrature
                                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                                end

                                if strcmp(pln.propDoseCalc.fineSampling.calcMode, 'fineSampling')
                                    % initialise empty dose array
                                    totalDose = zeros(size(currIx,1),1);

                                    if isfield(dij,'mLETDose')
                                        % calculate particle LET for bixel k on ray j of beam i
                                        depths = machine.data(energyIx).depths + machine.data(energyIx).offset;
                                        totalLET = zeros(size(currIx,1),1);
                                    end

                                    % run over components
                                    for c = 1:numOfSub(k)
                                        tmpDose = zeros(size(currIx,1),1);
                                        bixel = matRad_calcParticleDoseBixel(...
                                            currRadDepths(currIx(:,:,c),1,c), ...
                                            currRadialDist_sq(currIx(:,:,c),:,c), ...
                                            sigmaSub(k)^2, ...
                                            machine.data(energyIx));
                                        bixelDose = finalWeight(c,k).* bixel.physDose;

                                        tmpDose(currIx(:,:,c)) = bixelDose;
                                        totalDose = totalDose + tmpDose;

                                        if isfield(dij,'mLETDose')
                                            tmpLET = zeros(size(currIx,1),1);
                                            tmpLET(currIx(:,:,c)) = matRad_interp1(depths,machine.data(energyIx).LET,currRadDepths(currIx(:,:,c),1,c));
                                            totalLET = totalLET + tmpLET;
                                        end
                                    end

                                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix),1,totalDose,dij.doseGrid.numOfVoxels,1);
                                    if isfield(dij,'mLETDose')
                                        letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix),1,totalDose.*totalLET,dij.doseGrid.numOfVoxels,1);
                                    end
                                else
                                    % calculate particle dose for bixel k on ray j of beam i
                                    if isfield(pln,'propHeterogeneity') && pln.propHeterogeneity.calcHetero
                                        bixelDose = matRad_calcParticleDoseBixel(...
                                            currRadDepths(currIx), ...
                                            currRadialDist_sq(currIx), ...
                                            sigmaIni_sq, ...
                                            machine.data(energyIx), ...
                                            currHeteroCorrDepths, ...
                                            pln.propHeterogeneity, ...
                                            vTissueIndex_j(currIx));
                                    else
                                        bixelDose = matRad_calcParticleDoseBixel(...
                                            currRadDepths(currIx), ...
                                            currRadialDist_sq(currIx), ...
                                            sigmaIni_sq, ...
                                            machine.data(energyIx));
                                    end

                                    % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                                    %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                                    %Type               = 'dose';
                                    %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);

                                    % save dose for every bixel in cell array
                                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelDose.physDose,dij.doseGrid.numOfVoxels,1);

                                    if isfield(dij,'mLETDose')
                                        if isfield(bixelDose,'LET')
                                            bixelLET = bixelDose.LET;
                                        else
                                            % calculate particle LET for bixel k on ray j of beam i
                                            depths = machine.data(energyIx).depths + machine.data(energyIx).offset;
                                            bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,currRadDepths(currIx));
                                            bixelLET(isnan(bixelLET)) = 0;
                                        end

                                        % save LET for every bixel in cell array
                                        letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose.physDose,dij.doseGrid.numOfVoxels,1);
                                    end
                                end

                                % save alpha_p and beta_p radiosensititvy parameter for every bixel in cell array
                                if pln.bioParam.bioOpt

                                    if all(isfield(bixelDose,{'Z_Aij','Z_Bij'}))
                                        bixelAlphaDose =  bixelDose.L .* bixelDose.Z_Aij;
                                        bixelBetaDose  =  bixelDose.L .* bixelDose.Z_Bij;
                                    else
                                        if isfield(bixelDose,'LET') && pln.propHeterogeneity.modulateLET
                                            [bixelAlpha,bixelBeta] = pln.bioParam.calcLQParameter(currRadDepths,machine.data(energyIx),vTissueIndex_j(currIx,:),dij.ax(VdoseGrid(ix(currIx))),...
                                                dij.bx(VdoseGrid(ix(currIx))),dij.abx(VdoseGrid(ix(currIx))),bixelDose.LET);
                                        else
                                            [bixelAlpha,bixelBeta] = pln.bioParam.calcLQParameter(currRadDepths(currIx),machine.data(energyIx),vTissueIndex_j(currIx,:),...
                                                dij.ax(VdoseGrid(ix(currIx))),...
                                                dij.bx(VdoseGrid(ix(currIx))),...
                                                dij.abx(VdoseGrid(ix(currIx))));
                                        end
                                        bixelAlpha(isnan(bixelAlpha)) = 0;
                                        bixelBeta(isnan(bixelBeta)) = 0;

                                        bixelAlphaDose =  bixelDose.physDose .* bixelAlpha;
                                        bixelBetaDose  =  bixelDose.physDose .* sqrt(bixelBeta);
                                    end

                                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelAlphaDose,dij.doseGrid.numOfVoxels,1);
                                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen}  = sparse(VdoseGrid(ix(currIx)),1,bixelBetaDose,dij.doseGrid.numOfVoxels,1);

                                end
                            end
                        end
                    end

                    matRad_calcDoseFillDij;

                end % end bixels per ray

            end

        end %  end ray loop

    end % end beam loop


    % undo manipulation of isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(shiftScen,:);
    end

end % end shift scenario loop

dij = matRad_cleanDijScenarios(dij,pln,cst);

% Close Waitbar
if exist('figureWait') && ishandle(figureWait)
    delete(figureWait);
end