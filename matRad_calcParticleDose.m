function dij = matRad_calcParticleDose(ct,stf,pln,cst,param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper
% 
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   param:          (optional) structure defining additional parameter
%                   param.calcDoseDirect boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if exist('param','var')
    if ~isfield(param,'logLevel')
       param.logLevel = 1;
    end
    % default: dose influence matrix computation
   if ~isfield(param,'calcDoseDirect')
      param.calcDoseDirect = false;
   end
else
   param.calcDoseDirect = false;
   param.subIx          = [];
   param.logLevel       = 1;
end


if param.logLevel == 1
   % initialize waitbar
   figureWait = waitbar(0,'calculate dose influence matrix for particles...');
   % prevent closure of waitbar and show busy state
   set(figureWait,'pointer','watch');
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln, param);

% meta information for dij
dij.numOfBeams         = numel(stf);
dij.numOfVoxels        = prod(ct.cubeDim);
dij.resolution         = ct.resolution;
dij.dimensions         = ct.cubeDim;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);

if param.calcDoseDirect 
    numOfColumnsDij           = length(stf);
    numOfBixelsContainer = 1;
else
    numOfColumnsDij           = dij.totalNumOfBixels;
    numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
end

% set up arrays for book keeping
dij.bixelNum = NaN*ones(numOfColumnsDij,1);
dij.rayNum   = NaN*ones(numOfColumnsDij,1);
dij.beamNum  = NaN*ones(numOfColumnsDij,1);

% Allocate space for dij.physicalDose sparse matrix
for ctScen = 1:pln.multScen.numOfCtScen
    for ShiftScen = 1:pln.multScen.totNumShiftScen
        for RangeShiftScen = 1:pln.multScen.totNumRangeScen  
            
            if pln.multScen.scenMask(ctScen,ShiftScen,RangeShiftScen)
                dij.physicalDose{ctScen,ShiftScen,RangeShiftScen} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
            end
            
        end
    end
end


% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
doseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);


% if biological optimization considering a variable RBE is true then create alphaDose and betaDose containers and sparse matrices 
if pln.bioParam.bioOpt

    alphaDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.totNumShiftScen,pln.multScen.totNumRangeScen);
    
    for ctScen = 1:pln.multScen.numOfCtScen
        for ShiftScen = 1:pln.multScen.totNumShiftScen
            for RangeShiftScen = 1:pln.multScen.totNumRangeScen  
            
                if pln.multScen.scenMask(ctScen,ShiftScen,RangeShiftScen)
                    dij.mAlphaDose{ctScen,ShiftScen,RangeShiftScen}        = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
                    dij.mSqrtBetaDose{ctScen,ShiftScen,RangeShiftScen}     = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
                end
                
            end

        end

    end
    
end

% Only take voxels inside patient.
if isfield(param,'subIx') && ~isempty(param.subIx)
   V = param.subIx; 
else
   V = [cst{:,4}];
   V = unique(vertcat(V{:}));
end

% ignore densities outside of contours
eraseCtDensMask = ones(dij.numOfVoxels,1);
eraseCtDensMask(V) = 0;
for i = 1:dij.numOfScenarios
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% Convert CT subscripts to linear indices.
[yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
catch
   matRad_dispToConsole(['Could not find the following machine file: ' fileName ],param,'error'); 
end


if isfield(pln,'propDoseCalc') && ...
   isfield(pln.propDoseCalc,'calcLET') && ...
   pln.propDoseCalc.calcLET

  if isfield(machine.data,'LET')

    letDoseTmpContainer = cell(numOfBixelsContainer,pln.multScen.numOfCtScen,pln.multScen.numOfShiftScen,pln.multScen.numOfRangeShiftScen);
   
    for ctScen = 1:pln.multScen.numOfCtScen
        for ShiftScen = 1:pln.multScen.numOfShiftScen
            for RangeShiftScen = 1:pln.multScen.numOfRangeShiftScen  
            
                if pln.multScen.scenMask(ctScen,ShiftScen,RangeShiftScen)
                     dij.mLETDose{ctScen,ShiftScen,RangeShiftScen} = spalloc(prod(ct.cubeDim),numOfColumnsDij,1);
                end
                
            end
        end

    end
    
  else
    matRad_dispToConsole('LET not available in the machine data. LET will not be calculated.',param,'warning');
  end
end

% book keeping - this is necessary since pln is not used in optimization or
% matRad_calcCubes
if strcmp(pln.bioParam.model,'constRBE')
   dij.RBE = pln.bioParam.RBE;
end

% generates tissue class matrix for biological treatment planning and alpha_x, beta_x, vectors 
if pln.bioParam.bioOpt
   
    dij.alphaX = zeros(dij.numOfVoxels,1);
    dij.betaX  = zeros(dij.numOfVoxels,1);
    dij.abX    = zeros(dij.numOfVoxels,1);
    
    % show warning - it is important to use the same overlap priorities
    % during optimization. Otherwise biological dose calculation might be
    % incorrect when using different alpha/beta ratios throughout the patient
    matRad_dispToConsole(['matRad: Please use same priorities and tissue classes for optimization \n'],param,'info');
    %set overlap priorities
    cst = matRad_setOverlapPriorities(cst);

    % create radiosensitivity vectors
    for i = 1:size(cst,1)
        % check if cst is compatible 
        if ~isfield(cst{i,5},'alphaX') || ~isfield(cst{i,5},'betaX') 
           cst{i,5}.alphaX = 0.1;
           cst{i,5}.betaX  = 0.05;
           matRad_dispToConsole(['matRad: using default alpha_x and beta_x parameters for ' cst{i,2} ' \n'],param,'warning');
        end
        
        if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
            dij.alphaX(cst{i,4}{1}) = cst{i,5}.alphaX;
            dij.betaX(cst{i,4}{1})  = cst{i,5}.betaX;               
        end
  
    end
    
    % create alpha_x beta_x ratio vector
    dij.abX(dij.betaX>0) = dij.alphaX(dij.betaX>0)./dij.betaX(dij.betaX>0);
    
    % generates tissue class matrix for biological optimization
    vTissueIndex = zeros(size(V,1),1);
    
   if strcmp(pln.radiationMode,'carbon')

       for i = 1:size(cst,1)
           % find indices of structures related to V
           [~, row] = ismember(vertcat(cst{i,4}{:}),V,'rows'); 
           % check if base data contains alphaX and betaX
           if   isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
               % check if cst is compatiable 
               if ~isempty(cst{i,5})

                   IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                                    ismember(machine.data(1).betaX,cst{i,5}.betaX));

                   % check consitency of biological baseData and cst settings
                   if ~isempty(IdxTissue)
                       vTissueIndex(row) = IdxTissue;
                   else
                       matRad_dispToConsole('biological base data and cst inconsistent \n',param,'error');
                   end
               else
                   vTissueIndex(row) = 1;
                   matRad_dispToConsole(['matRad: tissue type of ' cst{i,2} ' was set to 1  \n'],param,'info');
               end
           else
               matRad_dispToConsole('base data is incomplement - alphaX and/or betaX is missing',param,'error');
           end

       end
       
       matRad_dispToConsole('done. \n',param,'info');
       
   end

end

ctScen  = 1;        % current ct scenario
matRad_dispToConsole('matRad: Particle dose calculation... \n',param,'info');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over all shift scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ShiftScen = 1:pln.multScen.totNumShiftScen
   
    % manipulate isocenter
    for k = 1:length(stf)
        stf(k).isoCenter = stf(k).isoCenter + pln.multScen.isoShift(ShiftScen,:);
    end
    
    matRad_dispToConsole(['shift scenario ' num2str(ShiftScen) ' of ' num2str(pln.multScen.totNumShiftScen) ':  \n'],param,'info');
    matRad_dispToConsole('matRad: Particle dose calculation... \n',param,'info');
    
    counter = 0;
    
    % compute SSDs
    stf = matRad_computeSSD(stf,ct,ctScen,param);

   for i = 1:numel(stf) % loop over all beams

       matRad_dispToConsole(['Beam ' num2str(i) ' of ' num2str(dij.numOfBeams) ':  \n'],param,'info');

       % remember beam and bixel number
       if param.calcDoseDirect
         dij.beamNum(i)    = i;
         dij.rayNum(i)     = i;
         dij.bixelNum(i)   = i;
       end

       bixelsPerBeam = 0;

       % convert voxel indices to real coordinates using iso center of beam i
       xCoordsV = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
       yCoordsV = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
       zCoordsV = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
       coordsV  = [xCoordsV yCoordsV zCoordsV];

       % Get Rotation Matrix
       % Do not transpose matrix since we usage of row vectors &
       % transformation of the coordinate system need double transpose

       % rotation around Z axis (gantry)
       rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

       % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
       rot_coordsV = coordsV*rotMat_system_T;

       rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
       rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
       rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

       % Calcualte radiological depth cube
       lateralCutoffRayTracing = 50;
       matRad_dispToConsole('matRad: calculate radiological depth cube...',param,'info');
       radDepthV = matRad_rayTracing(stf(i),ct,V,rot_coordsV,lateralCutoffRayTracing);
       matRad_dispToConsole('done. \n',param,'info');

       % get indices of voxels where ray tracing results are available
       radDepthIx = find(~isnan(radDepthV{1}));

       % limit rotated coordinates to positions where ray tracing is availabe
       rot_coordsV = rot_coordsV(radDepthIx,:);

       % Determine lateral cutoff
       matRad_dispToConsole('matRad: calculate lateral cutoff...',param,'info');
       cutOffLevel          = 0.99;
       visBoolLateralCutOff = 0;
       machine              = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),ctScen,visBoolLateralCutOff);
       matRad_dispToConsole('done. \n',param,'info');    

       for j = 1:stf(i).numOfRays % loop over all rays

           if ~isempty(stf(i).ray(j).energy)

               % find index of maximum used energy (round to keV for numerical reasons
               energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);

               maxLateralCutoffDoseCalc = max(machine.data(energyIx).LatCutOff.CutOff);

               % Ray tracing for beam i and ray j
               [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsV, ...
                                                        stf(i).sourcePoint_bev, ...
                                                        stf(i).ray(j).targetPoint_bev, ...
                                                        machine.meta.SAD, ...
                                                        radDepthIx, ...
                                                        maxLateralCutoffDoseCalc);
 

               % just use tissue classes of voxels found by ray tracer
               if pln.bioParam.bioOpt
                       vTissueIndex_j = vTissueIndex(ix,:);
               end

               for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                   counter       = counter + 1;
                   bixelsPerBeam = bixelsPerBeam + 1;

                   if param.logLevel <= 2
                      % Display progress and update text only 200 times
                      if mod(bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                              matRad_progress(bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                              floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                      end
                      if param.logLevel == 1
                         % update waitbar only 100 times if it is not closed
                         if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(figureWait)
                             waitbar(counter/dij.totalNumOfBixels,figureWait);
                         end
                      end
                   end
                   
                   % remember beam and bixel number
                   if ~param.calcDoseDirect
                      dij.beamNum(counter)  = i;
                      dij.rayNum(counter)   = j;
                      dij.bixelNum(counter) = k;
                   end

                   % find energy index in base data
                   energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
 
                   % create offset vector to account for additional offsets modelled in the base data and a potential 
                   % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                   % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow  
                   % interpolations in matRad_calcParticleDoseBixel() and avoid extrapolations.
                   offsetRadDepth = machine.data(energyIx).offset - stf(i).ray(j).rangeShifter(k).eqThickness;

                
                   for ctScen = 1:pln.multScen.numOfCtScen
                       for RangeShiftScen = 1:pln.multScen.totNumRangeScen 
                          
                          if pln.multScen.scenMask(ctScen,ShiftScen,RangeShiftScen)
                             
                            radDepths = radDepthV{ctScen}(ix);    
   
                            % manipulate radDepthCube for range scenarios 
                            if pln.multScen.relRangeShift(RangeShiftScen) ~= 0 || pln.multScen.absRangeShift(RangeShiftScen) ~= 0
                                    radDepths = radDepths +...                                                        % original cube
                                                radDepthV{ctScen}(ix)*pln.multScen.relRangeShift(RangeShiftScen) +... % rel range shift
                                                pln.multScen.absRangeShift(RangeShiftScen);                           % absolute range shift
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
                                matRad_dispToConsole('cutoff must be a value between 0 and 1',param,'error')
                            end

                            % empty bixels may happen during recalculation of error
                            % scenarios -> skip to next bixel
                            if ~any(currIx)
                                continue;
                            end

                            % adjust radDepth according to range shifter
                            currRadDepths = radDepths(currIx) + stf(i).ray(j).rangeShifter(k).eqThickness;

                            % calculate initial focus sigma
                            sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist(stf(i).ray(j).focusIx(k),:)', ...
                                                      machine.data(energyIx).initFocus.sigma(stf(i).ray(j).focusIx(k),:)',stf(i).ray(j).SSD{ctScen});
                            sigmaIni_sq = sigmaIni^2;

                            % consider range shifter for protons if applicable
                            if stf(i).ray(j).rangeShifter(k).eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                                % compute!
                                sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy, ...
                                                                   stf(i).ray(j).rangeShifter(k), ...
                                                                   stf(i).ray(j).SSD{ctScen});

                                % add to initial sigma in quadrature
                                sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                            end

                            % calculate particle dose for bixel k on ray j of beam i
                            bixelDose = matRad_calcParticleDoseBixel(...
                                currRadDepths, ...
                                radialDist_sq(currIx), ...
                                sigmaIni_sq, ...
                                machine.data(energyIx));       


                            % dij sampling is exluded for particles until we investigated the influence of voxel sampling for particles
                            %relDoseThreshold   =  0.02;   % sample dose values beyond the relative dose
                            %Type               = 'dose';
                            %[currIx,bixelDose] = matRad_DijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);

                            % save dose for every bixel in cell array
                            doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelDose,dij.numOfVoxels,1); 

                            if isfield(dij,'mLETDose')
                              % calculate particle LET for bixel k on ray j of beam i
                              depths   = machine.data(energyIx).depths + machine.data(energyIx).offset; 
                              bixelLET = matRad_interp1(depths,machine.data(energyIx).LET,radDepths(currIx)); 
                              bixelLET(isnan(bixelLET)) = 0;
                              % save LET for every bixel in cell array
                              letDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelLET.*bixelDose,dij.numOfVoxels,1);
                            end 

                            % save alpha_p and beta_p radiosensititvy parameter for every bixel in cell array 
                            if pln.bioParam.bioOpt

                               [bixelAlpha,bixelBeta] = pln.bioParam.calcLQParameter(radDepths(currIx),machine.data(energyIx),vTissueIndex_j(currIx,:),dij.alphaX(V(ix(currIx))),...
                                                                                                         dij.betaX(V(ix(currIx))),...
                                                                                                         dij.abX(V(ix(currIx))));
                                  
                                alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen} = sparse(V(ix(currIx)),1,bixelAlpha.*bixelDose,dij.numOfVoxels,1);
                                betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen}  = sparse(V(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.numOfVoxels,1);
                       
                            end
                         
                          end
                          
                       end

                   end
                   
                   % save computation time and memory by sequentially filling the
                   % sparse matrix dose.dij from the cell array
                   if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                      
                       for ctScen = 1:pln.multScen.numOfCtScen
                            for RangeShiftScen = 1:pln.multScen.totNumRangeScen
                                if ~any(currIx)
                                    continue;
                                end
                                if pln.multScen.scenMask(ctScen,ShiftScen,RangeShiftScen)
                  
                                      if param.calcDoseDirect
                                          if isfield(stf(1).ray(1),'weight') && numel(stf(i).ray(j).weight) >= k

                                              % score physical dose
                                              dij.physicalDose{ctScen,ShiftScen,RangeShiftScen}(:,i) = dij.physicalDose{ctScen,ShiftScen,RangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * doseTmpContainer{1,ctScen,ShiftScen,RangeShiftScen};

                                              if isfield(dij,'mLETDose') && pln.sampling
                                                   % score LETxDose matrices
                                                   dij.mLETDose{ctScen,ShiftScen,RangeShiftScen}(:,i) = dij.mLETDose{ctScen,ShiftScen,RangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * letDoseTmpContainer{1,ctScen,ShiftScen,RangeShiftScen}; 
                                              end
                                              
                                              if pln.bioParam.bioOpt
                                                   % score alphaxDose and sqrt(beta)xDose matrices
                                                   dij.mAlphaDose{ctScen,ShiftScen,RangeShiftScen}(:,i)    = dij.mAlphaDose{ctScen,ShiftScen,RangeShiftScen}(:,i)    + stf(i).ray(j).weight(k) * alphaDoseTmpContainer{1,ctScen,ShiftScen,RangeShiftScen};
                                                   dij.mSqrtBetaDose{ctScen,ShiftScen,RangeShiftScen}(:,i) = dij.mSqrtBetaDose{ctScen,ShiftScen,RangeShiftScen}(:,i) + stf(i).ray(j).weight(k) * betaDoseTmpContainer{1,ctScen,ShiftScen,RangeShiftScen};
                                              end
                                          else
                                              matRad_dispToConsole(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)],param,'error');
                                          end
                                      else

                                           % fill entire dose influence matrix
                                           dij.physicalDose{ctScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen}];

                                           if isfield(dij,'mLETDose')
                                               % fill entire LETxDose influence matrix
                                               dij.mLETDose{ctScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [letDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen}]; 
                                           end

                                           if pln.bioParam.bioOpt
                                               % fill entire alphaxDose influence and sqrt(beta)xDose influence matrices
                                               dij.mAlphaDose{ctScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter)    = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen}];
                                               dij.mSqrtBetaDose{ctScen,ShiftScen,RangeShiftScen}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,ctScen,ShiftScen,RangeShiftScen}];
                                           end

                                      end
                                end
                            end
                       end
                   end

               end % end bixels per ray

           end

       end % end ray loop
       
   end % end beam loop
   
   % manipulate isocenter
   for k = 1:length(stf)
       stf(k).isoCenter = stf(k).isoCenter - pln.multScen.isoShift(ShiftScen,:);
   end 
   
end % end shift scenario loop
        


try
  % wait 0.1s for closing all waitbars
  allWaitBarFigures = findall(0,'type','figure','tag','TMWWaitbar'); 
  delete(allWaitBarFigures);
  pause(0.1); 
catch
end
