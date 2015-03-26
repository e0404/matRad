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
dij.mBeta = spalloc(numel(ct),dij.totalNumOfBixels,1);

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
    mTissueClass = zeros(size(V,1),2);
    mTissueClass(:,1) = V;
    for i=1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(cst{i,8},V,'rows');        
        mTissueClass(row,2)=cst{i,9}.TissueClass;
    end
    
    SourceOfBioData = 'CNAO';%'GSI';%'CNAO';
    MultiClass = true;
    Counter = 0;
    switch SourceOfBioData
        % use existing four alpha curves for chordoma cells measured at the
        % GSI in Darmstadt
        case 'GSI'
            load('GSI_Chardoma_Carbon_BioData.mat');
            load('GSI_Chardoma_Carbon_BioData2.mat'); 
            
            if MultiClass == true
                EnergyBaseData = [baseData(:).energy];
                totalNumberOfEvaluations=length(BioData)*numel(EnergyBaseData);
                for currTissClass = 1:length(BioData)

                    for i=1:numel(EnergyBaseData)
                        [~, Index] = min(abs(BioData(currTissClass).energy-EnergyBaseData(i)));
                        vDepth = BioData(currTissClass).depths(:,Index);
                        for j = 1 : numel(vDepth)
                            dummyAlpha = zeros(numel(BioData(currTissClass).energy),1);
                            dummyBeta = zeros(numel(BioData(currTissClass).energy),1);
                            for k = 1 : numel(BioData(currTissClass).energy)
                                dummyAlpha(k) = interp1(BioData(currTissClass).depths(:,k), BioData(currTissClass).alpha(:,k), vDepth(j),'linear');
                                dummyBeta(k) = interp1(BioData(currTissClass).depths(:,k), BioData(currTissClass).beta(:,k), vDepth(j),'linear');
                            end
                            vAlpha(j) = interp1(BioData(currTissClass).energy, dummyAlpha, EnergyBaseData(i),'linear');
                            vBeta(j) = interp1(BioData(currTissClass).energy, dummyBeta, EnergyBaseData(i),'linear');
                        end

                        baseData(i).alpha(:,currTissClass) = vAlpha';
                        baseData(i).beta(:,currTissClass) = vBeta'; 
                        baseData(i).res_range(:,currTissClass) = vDepth;
                        Counter = Counter+1;
                        matRad_progress(Counter, totalNumberOfEvaluations);
                    end
                end
                
                
            else
                
                
                    % works for now just with one tissue class
                    vEnergiesMeasured = [stBioData{1,1}(1,1).Energy stBioData{1,1}(1,2).Energy ...
                                  stBioData{1,1}(1,3).Energy stBioData{1,1}(1,4).Energy];

                    % extract experimental measured biological data          
                    for i = 1:length(stBioData{1,1})
                        vAlphaMeasured(:,i) = stBioData{1,1}(1,i).Alpha;
                        vDepthMeasured(:,i) = stBioData{1,1}(1,i).Depths;
                    end

                    % extract available beam energies from baseData
                    for i = 1:numel(baseData)
                        vEnergies(i)=baseData(1,i).energy;
                    end
                    vEnergies = sort(vEnergies);

                    mBeta =zeros(size(vAlphaMeasured,1),1);
                    mAlpha=zeros(size(vAlphaMeasured,1),1);

                    for i=1:numel(vEnergies)
                        [~, Index] = min(abs(vEnergiesMeasured-vEnergies(i)));
                        vDepth = vDepthMeasured(:,Index);
                        for j = 1 : numel(vDepth)
                            dummyAlpha = zeros(numel(vEnergiesMeasured),1);
                            for k = 1 : numel(vEnergiesMeasured)
                                dummyAlpha(k) = interp1(vDepthMeasured(:,k), vAlphaMeasured(:,k), vDepth(j),'linear');
                            end
                            mAlpha(j) = interp1(vEnergiesMeasured, dummyAlpha, vEnergies(i),'linear');
                            mBeta(j) = 0.05;
                        end
                        baseData(i).res_range = vDepth;
                        baseData(i).alpha = mAlpha;
                        baseData(i).beta = mBeta;       
                    end

            end
            
        case 'CNAO'
            baseDataBio =matRadParseBioData([pwd filesep 'database_AB2']);
            TissueName = {'NormalTissue','chordoma','CHO','brainstem'};
            alpha_x = [0.3 0.1 0.18 0.053];
            beta_x = [0.1 0.05 0.028 0.0266];
             alpha_x = [0.1 0.1 0.1 0.1];
            beta_x = [0.05 0.05 0.05 0.05];
            
            if MultiClass == true
                 % assume that we have 3 tissue classes
                 NumTissueClasses =3;
                 for currTissClass = 1:NumTissueClasses
                     for j= 1:length(baseDataBio)
                            [~, index] = min(abs([baseData.energy]-baseDataBio(j).energy));
                            
                            PaddingValueAlpha = min(baseDataBio(j).dEdxA./baseData(index).Z);
                            tmp=interp1(baseDataBio(j).depths,baseDataBio(j).dEdxA,...
                                baseData(j).depths./10,'linear','extrap');
                            baseData(j).alpha(:,currTissClass) = tmp./baseData(index).Z;
                            
                            
                            tmp=interp1(baseDataBio(j).depths,baseDataBio(j).dEdxB,...
                                baseData(j).depths./10,'linear');
                            tmpNaN = (tmp./baseData(index).Z).^2;
                            a = find(~isnan(tmpNaN));
                            tmpNaN(isnan(tmpNaN))=tmpNaN(a(end));
                            baseData(j).beta(:,currTissClass)=tmpNaN;
                            
                            baseData(j).res_range(:,currTissClass) = (baseData(j).range - baseData(j).depths)./10;

                            
                            Counter = Counter+1;
                       
                             baseData(j).Tissue(currTissClass).Name = TissueName{1,currTissClass};
                             baseData(j).Tissue(currTissClass).Class = currTissClass;
                             baseData(j).Tissue(currTissClass).alphaX = alpha_x(currTissClass);
                             baseData(j).Tissue(currTissClass).betaX = beta_x(currTissClass);

                            
                            
                            matRad_progress(Counter, NumTissueClasses*length(baseDataBio));
                     end
                 end

            end
      
    end

      
%               figure,      subplot(221),plot(baseData(1).depths(:,1),baseData(1).beta(:,1),'LineWidth',3),title('115MeV','FontSize',20),grid on, hold on
%                                        %plot([min(baseData(1).depths(:,1)) max(baseData(1).depths(:,1))],[0.05 0.05],'color','r','LineWidth',3),hold off
%                            subplot(222),plot(baseData(20).depths(:,1),baseData(21).beta(:,1),'LineWidth',3),title('178MeV','FontSize',20), grid on, hold on
%                                        %plot([min(baseData(21).depths(:,1)) max(baseData(21).depths(:,1))],[0.05 0.05],'color','r','LineWidth',3),hold off
%                            subplot(223),plot(baseData(60).depths(:,1),baseData(60).beta(:,1),'LineWidth',3),title('277MeV','FontSize',20), grid on, hold on
%                                        %plot([min(baseData(60).depths(:,1)) max(baseData(60).depths(:,1))],[0.05 0.05],'color','r','LineWidth',3),hold off
%                            subplot(224),plot(baseData(121).depths(:,1),baseData(121).beta(:,1),'LineWidth',3),title('398MeV','FontSize',20), grid on, hold on
%                                        %plot([min(baseData(121).depths(:,1)) max(baseData(121).depths(:,1))],[0.05 0.05],'color','r','LineWidth',3),hold off

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
                    mTissueClass_j= mTissueClass(ix,:);
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
                    [bixelAlpha, bixelBeta] = mLQParams(...
                        radDepths(currIx),...
                        baseData(energyIx),...
                        mTissueClass_j(currIx,:),...
                        baseData);
                
                    
                end
   
                % Save dose for every bixel in cell array
                doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelDose,numel(ct),1);
                if pln.bioOptimization == true
                    alphaTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelAlpha,numel(ct),1);
                    betaTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} = sparse(V(ix(currIx)),1,bixelBeta,numel(ct),1);
                end
                % save computation time and memory by sequentially filling the 
                % sparse matrix dose.dij from the cell array
                if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
                    dij.dose(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    
                    if pln.bioOptimization == true
                        dij.mAlpha(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [alphaTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                        dij.mBeta(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
                    end
                end

            end
            
        end
        
    end
end
