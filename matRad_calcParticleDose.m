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

matRad_cfg =  MatRad_Config.instance();

% initialize waitbar
figureWait = waitbar(0,'calculate dose influence matrix for particles...');
% prevent closure of waitbar and show busy state
set(figureWait,'pointer','watch');


matRad_cfg.dispInfo('matRad: Particle dose calculation... \n');

if ~isfield(pln,'propHeterogeneity')
    pln.propHeterogeneity.calcHetero = false;
else
    if ~isfield(pln.propHeterogeneity,'calcHetero')
        pln.propHeterogeneity.calcHetero = matRad_cfg.propHeterogeneity.defaultCalcHetero;
    end
    if ~isfield(pln.propHeterogeneity,'useOriginalDepths')
        pln.propHeterogeneity.useOriginalDepths = matRad_cfg.propHeterogeneity.defaultUseOriginalDepths;
    end
    if ~isfield(pln.propHeterogeneity,'type')
        pln.propHeterogeneity.type = matRad_cfg.propHeterogeneity.defaultType;
    end
    if ~isfield(pln.propHeterogeneity,'modulateBioDose')
        pln.propHeterogeneity.modulateBioDose = matRad_cfg.propHeterogeneity.defaultModulateBioDose;
    end
end

if pln.propHeterogeneity.calcHetero
   cstOriginal = cst; 
end

% init dose calc
matRad_calcDoseInit;
    
% initialize lung heterogeneity correction and turn off if necessary files are missing
if pln.propHeterogeneity.calcHetero
    matRad_cfg.dispInfo('Heterogeneity correction enabled. \n');
    heteroCST = false;
    for i = 1:length(cst(:,1)) % scan cst for segmentation flagged for correction
        if isfield(cst{i,5},'HeterogeneityCorrection')
            heteroCST = true;
            break
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
if pln.propHeterogeneity.calcHetero
    
    lungVoxel = unique(cell2mat([cstOriginal{contains(cst(:,2),'lung','IgnoreCase',true),4}]'),'rows'); % get all lung voxel indices
    calcHeteroCorrStruct.cubeDim = ct.cubeDim;
    calcHeteroCorrStruct.numOfCtScen = pln.multScen.numOfCtScen;
    calcHeteroCorrStruct.resolution = ct.resolution;
    
    calcHeteroCorrStruct.cube = cell(1,ctScen);
    calcHeteroCorrStruct.cube(1,:) = {zeros(ct.cubeDim)};
    
    for shiftScen = 1:pln.multScen.numOfCtScen
        calcHeteroCorrStruct.cube{shiftScen}(lungVoxel(:,shiftScen)) = ct.cube{shiftScen}(lungVoxel(:,shiftScen));
    end
    
    
    if pln.propHeterogeneity.useOriginalDepths
        machine.data = matRad_checkBaseData(machine.data);
        if isstruct(machine.data(1).Z) && isfield(machine.data(1).Z,'profileORG')
            for i = 1:length(machine.data)
                machine.data(i).Z = machine.data(i).Z.profileORG;
            end
        elseif isstruct(machine.data(1).Z) && ~isfield(machine.data(1).Z,'profileORG')
            matRad_cfg.dispWarning('No original depths available in base data. Nothing changed.');
        else
            matRad_cfg.dispWarning('Base data depths are already in the desired format.');
        end
    end
    
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

if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'calcLET') 
    pln.propDoseCalc.calcLET = matRad_cfg.propDoseCalc.defaultCalcLET;
end

if  pln.propDoseCalc.calcLET
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
        matRad_cfg.dispWarning('\tLET not available in the machine data. LET will not be calculated.');
    end
end

%Toggles correction of small difference of current SSD to distance used
%in generation of base data (e.g. phantom surface at isocenter)
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc, 'airOffsetCorrection') 
    pln.propDoseCalc.airOffsetCorrection = true;
    
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
    matRad_cfg.dispWarning('Machine data does not contain BAMStoIsoDist. Using default value of %f mm\n.',BAMStoIsoDist);
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
    if strcmp(pln.bioParam.model,'LEM')
        if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
            
            matRad_cfg.dispInfo('\tloading biological base data...');
            
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
            
            matRad_cfg.dispInfo(' done.\n');
            
        else
            matRad_cfg.dispError('base data is incomplement - alphaX and/or betaX is missing');
        end
        
    else
        % parametrized biological models are based on the LET
        if ~isfield(machine.data,'LET')
            matRad_cfg.dispError('base data is incomplement - LET is missing');
        end
    end %  end is LEM model  
end


% lateral cutoff for raytracing and geo calculations
if ~isfield(pln,'propDoseCalc') || ~isfield(pln.propDoseCalc,'geometricCutOff')
    effectiveLateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;
end

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
        
        matRad_calcDoseInitBeam;
        
        % Calculate radiological depth cube for heterogeneity correction
        if pln.propHeterogeneity.calcHetero
            matRad_cfg.dispInfo('matRad: calculate radiological depth cube for heterogeneity correction...');
            heteroCorrDepthV = matRad_rayTracing(stf(i),calcHeteroCorrStruct,VctGrid,rot_coordsV,effectiveLateralCutoff);
            % HETERO interpolate hetero depth cube to dose grid resolution
            heteroCorrDepthV = matRad_interpRadDepth...
                (ct,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,heteroCorrDepthV);
            matRad_cfg.dispInfo('done.\n');
        end
        
        % Determine lateral cutoff
        matRad_cfg.dispInfo('\tmatRad: calculate lateral cutoff...');
        cutOffLevel = pln.propDoseCalc.lateralCutOff;
        visBoolLateralCutOff = 0;
        machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
        matRad_cfg.dispInfo('done.\n');
        
        for j = 1:stf(i).numOfRays % loop over all rays
            
            if ~isempty(stf(i).ray(j).energy)
                
                % find index of maximum used energy (round to keV for numerical reasons
                
                energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);
                
                maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);
                
                % Ray tracing for beam i and ray j
                [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                    stf(i).sourcePoint_bev, ...
                    stf(i).ray(j).targetPoint_bev, ...
                    machine.meta.SAD, ...
                    find(~isnan(radDepthVdoseGrid{1})), ...
                    maxLateralCutoffDoseCalc);
                
                if pln.propHeterogeneity.calcHetero
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
                    
                    
                    for ctScen = 1:pln.multScen.numOfCtScen
                        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
                            
                            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                
                                radDepths = radDepthVdoseGrid{ctScen}(ix);
                                
                                % manipulate radDepthCube for range scenarios
                                if pln.multScen.relRangeShift(rangeShiftScen) ~= 0 || pln.multScen.absRangeShift(rangeShiftScen) ~= 0
                                    radDepths = radDepths +...                                                       % original cube
                                        radDepthVdoseGrid{ctScen}(ix)*pln.multScen.relRangeShift(rangeShiftScen) +... % rel range shift
                                        pln.multScen.absRangeShift(rangeShiftScen);                                   % absolute range shift
                                    radDepths(radDepths < 0) = 0;
                                end
                                
                                
                                
                                % find depth depended lateral cut off
                                if cutOffLevel >= 1
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth;
                                elseif cutOffLevel < 1 && cutOffLevel > 0
                                    % perform rough 2D clipping
                                    currIx = radDepths <= machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                        radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);
                                    
                                    % peform fine 2D clipping
                                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                            (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
                                    end
                                else
                                    matRad_cfg.dispError('Lateral Cut-Off must be a value between 0 and 1!')
                                end
                                
                                % empty bixels may happen during recalculation of error
                                % scenarios -> skip to next bixel
                                if ~any(currIx)
                                    doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(dij.doseGrid.numOfVoxels,1);
                                    continue;
                                end

                                % adjust radDepth according to range shifter
                                if  pln.propDoseCalc.airOffsetCorrection
                                    currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness + dR;
                                    
                                    %sanity check due to negative corrections
                                    currRadDepths(currRadDepths < 0) = 0;
                                else
                                    currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;
                                end
                                
                                if pln.propHeterogeneity.calcHetero
                                    currHeteroCorrDepths = heteroCorrDepths(currIx);
                                end
                                    
                                % calculate initial focus sigma
                                sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist(stf(i).ray(j).focusIx(k),:)', ...
                                    machine.data(energyIx).initFocus.sigma(stf(i).ray(j).focusIx(k),:)',stf(i).ray(j).SSD);
                                sigmaIni_sq = sigmaIni^2;
                                
                                % consider range shifter for protons if applicable
                                if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')
                                    
                                    % compute!
                                    sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx), ...
                                        stf(i).ray(j).rangeShifter(k), ...
                                        stf(i).ray(j).SSD);
                                    
                                    % add to initial sigma in quadrature
                                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;
                                    
                                end
                                
                                % calculate particle dose for bixel k on ray j of beam i
                                if pln.propHeterogeneity.calcHetero
                                    bixelDose = matRad_calcParticleDoseBixel(...
                                        currRadDepths, ...
                                        radialDist_sq(currIx), ...
                                        sigmaIni_sq, ...
                                        machine.data(energyIx), ...
                                        currHeteroCorrDepths, ...
                                        pln.propHeterogeneity.type, ...
                                        pln.propHeterogeneity.modulateBioDose, vTissueIndex_j(currIx));
                                else
                                    bixelDose = matRad_calcParticleDoseBixel(...
                                        currRadDepths, ...
                                        radialDist_sq(currIx), ...
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
                                    % calculate particle LET for bixel k on ray j of beam i
                                    depths   = machine.data(energyIx).depths + machine.data(energyIx).offset;
                                    bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,currRadDepths);
                                    bixelLET(isnan(bixelLET)) = 0;
                                    % save LET for every bixel in cell array
                                    letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,shiftScen,rangeShiftScen} = sparse(VdoseGrid(ix(currIx)),1,bixelLET.*bixelDose.physDose,dij.doseGrid.numOfVoxels,1);
                                end
                                
                                % save alpha_p and beta_p radiosensititvy parameter for every bixel in cell array
                                if pln.bioParam.bioOpt
                                    
                                    if all(isfield(bixelDose,{'Z_Aij','Z_Bij'}))
                                        bixelAlphaDose =  bixelDose.L .* bixelDose.Z_Aij;
                                        bixelBetaDose  =  bixelDose.L .* bixelDose.Z_Bij;
                                    else
                                    [bixelAlpha,bixelBeta] = pln.bioParam.calcLQParameter(currRadDepths,machine.data(energyIx),vTissueIndex_j(currIx,:),dij.ax(VdoseGrid(ix(currIx))),...
                                        dij.bx(VdoseGrid(ix(currIx))),dij.abx(VdoseGrid(ix(currIx))));
                                    
                                    bixelAlphaDose =  bixelDose.physDose .* bixelAlpha;
                                    bixelBetaDose  =  bixelDose.physDose .* sqrt(bixelBeta);
                                    end

                                    bixelAlpha(isnan(bixelAlpha)) = 0;
                                    bixelBeta(isnan(bixelBeta)) = 0;
                                    
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

%Close Waitbar
if ishandle(figureWait)
    delete(figureWait);
end