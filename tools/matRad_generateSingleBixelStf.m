function stf = matRad_generateSingleBixelStf(ct,cst,pln)
% 
% call
%   stf = matRad_generateSingleBixelStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating single bixel stf struct... ');

if nargin < 4
    visMode = 0;
end

if numel(pln.propStf.gantryAngles) ~= numel(pln.propStf.couchAngles)
    matRad_cfg.dispError('Inconsistent number of gantry and couch angles.');
end

if pln.propStf.bixelWidth < 0 || ~isfinite(pln.propStf.bixelWidth)
   matRad_cfg.dispError('bixel width (spot distance) needs to be a real number [mm] larger than zero.');
end

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine '.mat'];
try
   load([matRad_cfg.matRadRoot filesep 'basedata' filesep fileName]);
   SAD = machine.meta.SAD;
catch
   matRad_cfg.dispError('Could not find the following machine file: %s',fileName); 
end

if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
      
    availableEnergies = [machine.data.energy];
    availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
    
    if sum(availablePeakPos<0)>0
       matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file') 
    end
    %clear machine;
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

% Define steering file like struct. Prellocating for speed.
stf = struct;

% find set of energyies with adequate spacing
if ~isfield(pln.propStf, 'longitudinalSpotSpacing')
    longitudinalSpotSpacing = matRad_cfg.propStf.defaultLongitudinalSpotSpacing;
else
    longitudinalSpotSpacing = pln.propStf.longitudinalSpotSpacing;
end

if ~isfield(pln.propStf,'useRangeShifter')
    pln.propStf.useRangeShifter = false;
end

%Get Isocenter voxel as target
if ~any(isfield(ct,{'x','y','z'}))
    ct.x = ct.resolution.x*[1:ct.cubeDim(2)]-ct.resolution.x/2;
    ct.y = ct.resolution.y*[1:ct.cubeDim(1)]-ct.resolution.y/2;
    ct.z = ct.resolution.z*[1:ct.cubeDim(3)]-ct.resolution.z/2;
end

%xVox = ct.x + ct.resolution.x/2;
%yVox = ct.y + ct.resolution.y/2;
%zVox = ct.z + ct.resolution.z/2;

xVox = ct.resolution.x*[1:ct.cubeDim(2)]-ct.resolution.x/2;
yVox = ct.resolution.y*[1:ct.cubeDim(1)]-ct.resolution.y/2;
zVox = ct.resolution.z*[1:ct.cubeDim(3)]-ct.resolution.z/2;

% loop over all angles
for i = 1:length(pln.propStf.gantryAngles)   

    % Save meta information for treatment plan
    stf(i).gantryAngle          = pln.propStf.gantryAngles(i);
    stf(i).couchAngle           = pln.propStf.couchAngles(i);
    stf(i).bixelWidth           = pln.propStf.bixelWidth;
    stf(i).radiationMode        = pln.radiationMode;
    stf(i).SAD                  = SAD;
    stf(i).isoCenter            = pln.propStf.isoCenter(i,:);
    stf(i).numOfRays            = 1;
    stf(i).numOfBixelsPerRay    = 1;
    stf(i).totalNumOfBixels     = 1;
    
    x = floor(matRad_interp1(xVox,[1:ct.cubeDim(2)]',stf.isoCenter(1)));
    y = floor(matRad_interp1(yVox,[1:ct.cubeDim(1)]',stf.isoCenter(2)));
    z = floor(matRad_interp1(zVox,[1:ct.cubeDim(3)]',stf.isoCenter(3)));
    
    %Voxel index of Isocenter
    isoIx = [y x z];
    
    
    % generate voi cube for targets
    voiTarget    = zeros(ct.cubeDim);
    voiTarget(isoIx(1),isoIx(2),isoIx(3)) = 1;
    %adds = unique(perms([1 0 0]),'rows');
    %adds = [adds; -adds];
    %for p = 1:size(adds,1)
    %    ix = isoIx + adds(i,:);
    %    voiTarget(ix(1),ix(2),ix(3)) = 1;
    %end
    
    %voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution);
        
    % Get the (active) rotation matrix. We perform a passive/system 
    % rotation with row vector coordinates, which would introduce two 
    % inversions / transpositions of the matrix, thus no changes to the
    % rotation matrix are necessary
    rotMat_system_T = matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i));
    rotMat_vectors_T = transpose(rotMat_system_T);
    
    stf(i).sourcePoint_bev  = [0 -SAD 0];
    stf(i).sourcePoint      = stf(i).sourcePoint_bev*rotMat_vectors_T;  
    
    stf(i).ray.rayPos_bev = [0 0 0];
    stf(i).ray.targetPoint_bev = [0 SAD 0];
    
    stf(i).ray.rayPos = stf.isoCenter;
    stf(i).ray.targetPoint = [0 SAD 0] * rotMat_vectors_T;
    
    
    % find appropriate energies for particles
    if strcmp(stf(i).radiationMode,'protons') || strcmp(stf(i).radiationMode,'carbon')
        stf.longitudinalSpotSpacing = longitudinalSpotSpacing;
              
        % ray tracing necessary to determine depth of the target
        [alphas,l,rho,d12,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
                             ct.resolution, ...
                             stf(i).sourcePoint, ...
                             stf(i).ray.targetPoint, ...
                             [{ct.cube{1}} {voiTarget}]);
                         
        %Used for generic range-shifter placement                  
        ctEntryPoint = alphas(1) * d12;
        
        % target hit
        if sum(rho{2}) > 0
            
            % compute radiological depths
            % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
            radDepths = cumsum(l .* rho{1});
            
            % find target entry & exit
            diff_voi    = diff([rho{2}]);
            targetEntry = radDepths(diff_voi == 1);
            targetExit  = radDepths(diff_voi == -1);
            
            if numel(targetEntry) ~= numel(targetExit)
                matRad_cfg.dispError('Inconsistency during ray tracing. Please check correct assignment and overlap priorities of structure types OAR & TARGET.');
            end
            
            raShiThickness = 50;
            
            % Save energies in stf struct
            bestPeakPos = mean([targetExit,targetEntry]) + raShiThickness;
            [~,closest] = min(abs(availablePeakPos - bestPeakPos));
            stf(i).ray.energy = availableEnergies(closest);
            
            
            % book keeping & calculate focus index
            currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
                machine.meta.LUT_bxWidthminFWHM(2,:)',...
                pln.propStf.bixelWidth, ...
                machine.meta.LUT_bxWidthminFWHM(2,end));            
            [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                repmat(stf(i).ray.energy,length([machine.data]),1))));
            
            % get for each spot the focus index
            stf(i).ray.focusIx = find(machine.data(vEnergyIx).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
            
            if pln.propStf.useRangeShifter
            
                %Include range shifter data
                stf(i).ray.rangeShifter.ID = 1;
                stf(i).ray.rangeShifter.eqThickness = raShiThickness;
            
                %Place range shifter 2 times the range away from isocenter, but
                %at least 10 cm
                sourceRaShi = round(ctEntryPoint - 2*raShiThickness,-1); %place a little away from entry, round to cms to reduce number of unique settings;
                stf(i).ray.rangeShifter.sourceRashiDistance = sourceRaShi;
            else
                stf(i).ray.rangeShifter.ID = 0;
                stf(i).ray.rangeShifter.eqThickness = 0;
                stf(i).ray.rangeShifter.sourceRashiDistance = 0;
            end
            
            
        else % target not hit
            stf(i).ray               = [];
            stf(i).numOfBixelsPerray = [];
        end
        
    elseif strcmp(stf(i).radiationMode,'photons')
        
        % book keeping for photons
        stf(i).ray.energy = machine.data.energy;
        
    else
        matRad_cfg.dispError('Error generating stf struct: invalid radiation modality.');
    end
end    



