function dij = matRad_calcParticleDose(ct,stf,pln,cst,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad particle dose calculation wrapper
% 
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst,visBool)
%
% input
%   ct:         ct cube
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%   cst:        matRad cst struct
%   visBool:    toggle on/off visualization (optional)
%
% output
%   dij:        matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if visBool not set toogle off visualization
if nargin < 5
    visBool = 0;
end

% initialize waitbar
figureWait = waitbar(0,'calculate particle-ij matrice(s)...');

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = ct.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.physicalDose sparse matrix
dij.physicalDose = spalloc(numel(ct.cube),dij.totalNumOfBixels,1);

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% Allocate memory for dose_temp cell array
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
doseTmpContainer = cell(numOfBixelsContainer,1);
if pln.bioOptimization == true 
    alphaDoseTmpContainer = cell(numOfBixelsContainer,1);
    betaDoseTmpContainer  = cell(numOfBixelsContainer,1);
    dij.mAlphaDose        = spalloc(numel(ct.cube),dij.totalNumOfBixels,1);
    dij.mSqrtBetaDose     = spalloc(numel(ct.cube),dij.totalNumOfBixels,1);
end

% Only take voxels inside patient.
V = unique([cell2mat(cst(:,4))]);

% Convert CT subscripts to linear indices.
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct.cube),V);

xCoordsV = xCoordsV(:)*ct.resolution.x-pln.isoCenter(1);
yCoordsV = yCoordsV(:)*ct.resolution.y-pln.isoCenter(2);
zCoordsV = zCoordsV(:)*ct.resolution.z-pln.isoCenter(3);
coordsV  = [xCoordsV yCoordsV zCoordsV];

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]); 
end

% generates tissue class matrix for biological optimization
if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
        && strcmp(pln.radiationMode,'carbon')
    fprintf('matRad: loading biological base data... ');
    mTissueClass = zeros(size(V,1),1);
    for i = 1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(cst{i,4},V,'rows');  
        if ~isempty(cst{i,5}) && isfield(cst{i,5},'TissueClass')
            mTissueClass(row) = cst{i,5}.TissueClass;
        else
            mTissueClass(row) = 1;
            fprintf(['matRad: tissue type of ' cst{i,2} ' was set to 1 \n']);
        end
        
        % check consitency of biological baseData and cst settings
        baseDataAlphaBetaRatios =  reshape([machine.data(:).alphaBetaRatio],numel(machine.data(1).alphaBetaRatio),size(machine.data,2));
        if norm(baseDataAlphaBetaRatios(cst{i,5}.TissueClass,:) - cst{i,5}.alphaX/cst{i,5}.betaX)>0
            error('biological base data and cst inconsistent\n');
        end
        
    end
    fprintf('...done \n');
end

% source position in beam's eye view.
sourcePoint_bev = [0 -pln.SAD 0];

% determine lateral cutoff
fprintf('matRad: calculate lateral cutoff... ');
cutOffLevel = 0.99;
visBoolLateralCutOff = 0;
machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,visBoolLateralCutOff);
fprintf('...done \n');

fprintf('matRad: Particle dose calculation... ');
counter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    % Set gantry and couch rotation matrices according to IEC 61217
    % Use transpose matrices because we are working with row vectros
    
    % rotation around Z axis (gantry)
    inv_rotMx_XY_T = [ cosd(-pln.gantryAngles(i)) sind(-pln.gantryAngles(i)) 0;
                      -sind(-pln.gantryAngles(i)) cosd(-pln.gantryAngles(i)) 0;
                                                0                          0 1];
    
    % rotation around Y axis (couch)
    inv_rotMx_XZ_T = [cosd(-pln.couchAngles(i)) 0 -sind(-pln.couchAngles(i));
                                              0 1                         0;
                      sind(-pln.couchAngles(i)) 0  cosd(-pln.couchAngles(i))];
                  
    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV = coordsV*inv_rotMx_XZ_T*inv_rotMx_XY_T;
    
    rot_coordsV(:,1) = rot_coordsV(:,1)-sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-sourcePoint_bev(3);
    
    for j = 1:stf(i).numOfRays % loop over all rays
        
        if ~isempty(stf(i).ray(j).energy)
        
            % find index of maximum used energy (round to keV for numerical
            % reasons
            energyIx = max(round2(stf(i).ray(j).energy,4)) == round2([machine.data.energy],4);
            
            
            maxLateralCutoff = max(machine.data(energyIx).LatCutOff.CutOff);
            
            
            % Ray tracing for beam i and ray j
            [ix,radDepths,~,latDistsX,latDistsZ] = matRad_calcRadGeoDists(ct.cube, ...
                                                        V,...
                                                        pln.isoCenter, ...
                                                        rot_coordsV, ...
                                                        ct.resolution, ...
                                                        stf(i).sourcePoint, ...
                                                        stf(i).ray(j).targetPoint, ...
                                                        sourcePoint_bev,...
                                                        stf(i).ray(j).targetPoint_bev, ...
                                                        coordsV, ...
                                                        maxLateralCutoff, ...
                                                        visBool);
            
            radialDist_sq = latDistsX.^2 + latDistsZ.^2;    
            
            % just use tissue classes of voxels found by ray tracer
            if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                 && strcmp(pln.radiationMode,'carbon')
                    mTissueClass_j= mTissueClass(ix,:);
            end
            
            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                % Display progress
                matRad_progress(counter,dij.totalNumOfBixels);
                % update waitbar only 100 times
                if mod(counter,round(dij.totalNumOfBixels/100)) == 0
                    waitbar(counter/dij.totalNumOfBixels);
                end
                % remember beam and  bixel number
                dij.beamNum(counter)  = i;
                dij.rayNum(counter)   = j;
                dij.bixelNum(counter) = k;

                % find energy index in base data
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                
                % find depth depended lateral cut off
                if cutOffLevel >= 1
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset;
                elseif cutOffLevel < 1 && cutOffLevel > 0
                    % perform rough 2D clipping
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset & ...
                         radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                    % peform fine 2D clipping  
                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                        currIx(currIx) = interp1(machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset,...
                            machine.data(energyIx).LatCutOff.CutOff.^2, radDepths(currIx)) >= radialDist_sq(currIx);
                    end
                else
                    error('cutoff must be a value between 0 and 1')
                end
                 
                % calculate particle dose for bixel k on ray j of beam i
                bixelDose = matRad_calcParticleDoseBixel(...
                    radDepths(currIx),...
                    radialDist_sq(currIx),...
                    machine.data(energyIx)); 
                
                % Save dose for every bixel in cell array
                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelDose,numel(ct.cube),1);
                            
                if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                    && strcmp(pln.radiationMode,'carbon')
                    % calculate alpha and beta values for bixel k on ray j of                  
                    [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                        radDepths(currIx),...
                        mTissueClass_j(currIx,:),...
                        machine.data(energyIx));
                
                    alphaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelAlpha.*bixelDose,numel(ct.cube),1);
                    betaDoseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,sqrt(bixelBeta).*bixelDose,numel(ct.cube),1);
                end
                
                % save computation time and memory by sequentially filling the
                % sparse matrix dose.dij from the cell array
                if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                    dij.physicalDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    
                    if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                            && strcmp(pln.radiationMode,'carbon')
                        dij.mAlphaDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                        dij.mSqrtBetaDose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    end
                end

            end
            
        end
        
    end
end
close(figureWait);
