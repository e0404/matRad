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
%   -
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

% meta information for dij
dij.numOfBeams         = pln.numOfBeams;
dij.numOfVoxels        = pln.numOfVoxels;
dij.resolution         = pln.resolution;
dij.numOfRaysPerBeam   = [stf(:).numOfRays];
dij.totalNumOfRays     = sum(dij.numOfRaysPerBeam);
dij.totalNumOfBixels   = sum([stf(:).totalNumOfBixels]);
dij.dimensions         = pln.voxelDimensions;

% set up arrays for book keeping
dij.bixelNum = NaN*ones(dij.totalNumOfRays,1);
dij.rayNum   = NaN*ones(dij.totalNumOfRays,1);
dij.beamNum  = NaN*ones(dij.totalNumOfRays,1);

% Allocate space for dij.dose sparse matrix
dij.dose = spalloc(numel(ct),dij.totalNumOfBixels,1);
dij.mAlpha = spalloc(numel(ct),dij.totalNumOfBixels,1);
dij.mBeta = 0.05;

% Allocate memory for dose_temp cell array
numOfBixelsContainer = ceil(dij.totalNumOfBixels/10);
doseTmpContainer = cell(numOfBixelsContainer,1);
if pln.bioOptimization == true 
    alphaTmpContainer = cell(numOfBixelsContainer,1);
end
% Only take voxels inside patient.
V = unique([cell2mat(cst(:,8))]);

% Convert CT subscripts to linear indices.
[yCoordsV, xCoordsV, zCoordsV] = ind2sub(size(ct),V);

xCoordsV = (xCoordsV(:)-0.5)*pln.resolution(1)-pln.isoCenter(1);
yCoordsV = (yCoordsV(:)-0.5)*pln.resolution(2)-pln.isoCenter(2);
zCoordsV = (zCoordsV(:)-0.5)*pln.resolution(3)-pln.isoCenter(3);
coords_inside=[xCoordsV yCoordsV zCoordsV];


% load protonBaseData
if strcmp(pln.radiationMode,'protons')
    load protonBaseData;
elseif strcmp(pln.radiationMode,'carbon')
    load carbonBaseData;
end


% generates tissue class matrix for biological optimization
% and initializes alpha/beta interpolants
if pln.bioOptimization == true
    fprintf('matRad: Creating biological tissue interpolant... ');
    mT = zeros(size(V,1),2);
    mT(:,1) = V;
    for i=1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(cst{i,8},V,'rows');        
        mT(row,2)=cst{i,9}.TissueClass;
    end
    
    load('GSI_Chardoma_Carbon_BioData.mat');
    
    
    
    
    
    tTEnergies = [stBioData{1,1}(1,1).Energy stBioData{1,1}(1,2).Energy ...
                  stBioData{1,1}(1,3).Energy stBioData{1,1}(1,4).Energy];
    tTAlpha = zeros(81,4);
    tTAlpha(:,1)=stBioData{1,1}(1,1).Alpha;
    tTAlpha(:,2)=stBioData{1,1}(1,2).Alpha;
    tTAlpha(:,3)=stBioData{1,1}(1,3).Alpha;
    tTAlpha(:,4)=stBioData{1,1}(1,4).Alpha;
    tDepth = zeros(81,4);
    tDepth(:,1)=stBioData{1,1}(1,1).Depths;
    tDepth(:,2)=stBioData{1,1}(1,2).Depths;
    tDepth(:,3)=stBioData{1,1}(1,3).Depths;
    tDepth(:,4)=stBioData{1,1}(1,4).Depths;
    
    
    vEnergies = zeros(numel(baseData),1);
    for i = 1:numel(baseData)
        vEnergies(i)=baseData(1,i).energy;
    end
    vEnergies = sort(vEnergies);
    vDepth = linspace(min(tDepth(:)),max(tDepth(:)),300);
    mAlpha = zeros(numel(vDepth),numel(vEnergies));
    mDepth =zeros(numel(vDepth),numel(vEnergies));
   
    for i=1:numel(vEnergies)
        for IX = 1 : numel(vDepth)

            dummyAlpha = zeros(numel(tTEnergies),1);

            for JX = 1 : numel(tTEnergies)
                dummyAlpha(JX) = interp1(tDepth(:,JX), tTAlpha(:,JX), vDepth(IX),'linear');
            end

            mAlpha(IX,i) = interp1(tTEnergies, dummyAlpha, vEnergies(i),'linear');

        end
        mDepth(:,i)=vDepth;
    end

    BioInterp.tTEnergies =vEnergies;
    BioInterp.tTAlpha=mAlpha;
    BioInterp.tDepth=mDepth;
    
    
%      figure, subplot(221),plot(tDepth(:,1),tTAlpha(:,1)),title('88MeV'),
%             subplot(222),plot(tDepth(:,2),tTAlpha(:,2)),title('178MeV'),
%             subplot(223),plot(tDepth(:,3),tTAlpha(:,3)),title('276MeV'),
%             subplot(224),plot(tDepth(:,4),tTAlpha(:,4)),title('430MeV')
%     
%     figure,
%     for i=1:size(mDepth,2)
%         str = sprintf('Energy %d',vEnergies(i));
%         plot(mDepth(:,i),mAlpha(:,i)),title(str);
%         waitforbuttonpress
%     end


   
        fprintf('...done \n');
end


% It make a meshgrid with CT position in millimeter for calculate
% geometrical distances
[X_geo,Y_geo,Z_geo] = meshgrid(pln.resolution(1)*(0.5:1:size(ct,1)),...
    pln.resolution(2)*(0.5:1:size(ct,2)),pln.resolution(3)*(0.5:1:size(ct,3)));

% take only voxels inside patient
X_geo = X_geo(V);
Y_geo = Y_geo(V);
Z_geo = Z_geo(V);

% source position in beam's eye view.
sourcePoint_bev = [0 -pln.SAD 0];

counter = 0;

fprintf('matRad: Particle dose calculation... ');

if strcmp(pln.radiationMode,'protons')
    mLQParams = @(FreeParameter) matRad_ProtonLQParameter(FreeParameter,0);
elseif strcmp(pln.radiationMode,'carbon')
    mLQParams = @(vRadDepths,sEnergy,mT,Interp) matRad_CarbonLQParameter(vRadDepths,sEnergy,mT,Interp,0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:dij.numOfBeams; % loop over all beams
    
    %SET GANTRY AND COUCH ROTATION MATRICES ACCORDING IEC 60601 STANDARD FOR LINACS
    % Note: Signs for the following 2 matrices works for a fixed beam and
    % rotary CT.
    
    % Rotation around Z axis (gantry movement)
    rotMx_XY = [ cosd(pln.gantryAngles(i)) -sind(pln.gantryAngles(i)) 0;
        sind(pln.gantryAngles(i)) cosd(pln.gantryAngles(i)) 0;
        0 0 1];
    
    % Rotation around Y axis (Couch movement)
    rotMx_XZ = [cosd(pln.couchAngles(i)) 0 -sind(pln.couchAngles(i));
        0 1 0;
        sind(pln.couchAngles(i)) 0 cosd(pln.couchAngles(i))];
    
    % ROTATE VOI'S CT COORDINATES. First applies couch rotation and then
    % gantry. It is important to note matrix multiplication is not "commutative",
    % you cannot switch the order of the factors and expect to end up with the same result.
    
    % Rotate coordinates around Y axis (1st couch movement) and then Z axis
    % (2nd gantry movement)
    
    [rot_coords] = coords_inside*rotMx_XZ*rotMx_XY;
    
    
    for j = 1:stf(i).numOfRays % loop over all rays
        
        if ~isempty(stf(i).ray(j).energy)
        
            % set lateral cutoff for calculation of geometric distances
            lateralCutoff = 3*baseData(find(max(stf(i).ray(j).energy) == [baseData.energy])).sigma(end);

            % Ray tracing for beam i and ray j
            [ix,radDepths,~,latDistsX,latDistsZ] = matRad_calcRadGeoDists(ct,V,...
                    pln.isoCenter,rot_coords,pln.resolution,stf(i).sourcePoint,...
                    stf(i).ray(j).targetPoint,sourcePoint_bev,...
                    stf(i).ray(j).targetPoint_bev,X_geo,Y_geo,Z_geo,lateralCutoff,visBool);
            
            radialDist_sq = latDistsX.^2 + latDistsZ.^2;    
            
            % just use tissue classes of voxels found by ray tracer
            if pln.bioOptimization == true 
                    mT_j= mT(ix,:);
            end
            
            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;
                % Display progress
                matRad_progress(counter,dij.totalNumOfBixels);

                % remember beam and  bixel number
                dij.beamNum(counter)  = i;
                dij.rayNum(counter)   = j;
                dij.bixelNum(counter) = k;

                % find energy index in base data
                energyIx = find(stf(i).ray(j).energy(k) == [baseData.energy]);

                % find indices
                currIx = radDepths <= baseData(energyIx).depths(end) & ...
                         radialDist_sq <= 9*baseData(energyIx).sigma(end)^2;

                % calculate particle dose for bixel k on ray j of beam i
                bixelDose = matRad_calcParticleDoseBixel(...
                    radDepths(currIx),...
                    radialDist_sq(currIx),...
                    baseData(energyIx));
                
            
                if pln.bioOptimization == true 
                    % calculate alpha and beta values for bixel k on ray j of
                    % beam i - call duration 0.0020s                    
                    [bixelAlpha, ~] = mLQParams(...
                        radDepths(currIx),...
                        baseData(energyIx),...
                        mT_j(currIx,:),...
                        BioInterp);
                
                    
                end
   
                % Save dose for every bixel in cell array
                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelDose,numel(ct),1);
                if pln.bioOptimization == true
                    alphaTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelAlpha,numel(ct),1);
                end
                % save computation time and memory by sequentially filling the 
                % sparse matrix dose.dij from the cell array
                if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                    dij.dose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    
                    if pln.bioOptimization == true
                        dij.mAlpha(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [alphaTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    end
                end

            end
            
        end
        
    end
end
