function stf = matRad_StfVMATPost(stf,pln,masterRayPosBEV,masterTargetPointBEV,SAD,machine)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad steering information post-processing for VMAT
%
% call
%   stf = matRad_StfVMATPost(stf,pln)
%
% input
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%After all steering file information is completed, loop over
%initialized gantry angles.  All children and subchildren of these angles should
%have ray positions given by the union of their own ray positions and
%the ray positions of the parent transformed to their gantry angle.
%This is so that: (1) the target is still totally in the FOV of each
%angle; and (2) the parent can give segments to the children during
%initial segmentation and DAO.

fprintf('matRad: VMAT post-processing (1/2)... ');

for i = 1:length(pln.propStf.gantryAngles)
    
    % Determine which FMO beam the current beam belongs to
    [~,stf(i).propVMAT.beamParentFMOIndex] = min(abs(pln.propStf.FMOGantryAngles-pln.propStf.gantryAngles(i)));
    stf(i).propVMAT.beamParentGantryAngle = pln.propStf.FMOGantryAngles(stf(i).propVMAT.beamParentFMOIndex);
    stf(i).propVMAT.beamParentIndex = find(pln.propStf.gantryAngles == stf(i).propVMAT.beamParentGantryAngle);
    
    % Indicate if this beam is to be included in DOA/FMO or not. All beams 
    % are still considered in dose calc for objective function in DAO
    stf(i).propVMAT.FMOBeam = any(pln.propStf.FMOGantryAngles==pln.propStf.gantryAngles(i));
    stf(i).propVMAT.DAOBeam = any(pln.propStf.DAOGantryAngles==pln.propStf.gantryAngles(i));
    
    %% Determine different angle borders
    
    % doseAngleBorders are the angular borders over which dose is deposited
    if i == 1
        
        stf(i).propVMAT.doseAngleBorders = ([pln.propStf.gantryAngles(i) pln.propStf.gantryAngles(i+1)]+pln.propStf.gantryAngles(i))/2;
    elseif i == length(pln.propStf.gantryAngles)
        
        stf(i).propVMAT.doseAngleBorders = ([pln.propStf.gantryAngles(i-1) pln.propStf.gantryAngles(i)]+pln.propStf.gantryAngles(i))/2;
    else
        
        stf(i).propVMAT.doseAngleBorders = ([pln.propStf.gantryAngles(i-1) pln.propStf.gantryAngles(i+1)]+pln.propStf.gantryAngles(i))/2;
    end
    
    stf(i).propVMAT.doseAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.doseAngleBorders(1) stf(i).propVMAT.doseAngleBorders(2)-stf(i).gantryAngle];
    stf(i).propVMAT.doseAngleBordersDiff = sum(stf(i).propVMAT.doseAngleBorderCentreDiff);
    
    %Assign beam to its Parent, either as child (optimized) or subchild
    %(interpolated)
    if stf(i).propVMAT.DAOBeam
        if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles)
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = 0;
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles = nan(1000,1);
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex = nan(1000,1);
        end
        
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren+1;
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = pln.propStf.gantryAngles(i);
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamChildren) = i;
        
        %optAngleBorders are the angular borders over which an optimized control point
        %has influence
        DAOIndex = find(pln.propStf.DAOGantryAngles == pln.propStf.gantryAngles(i));
        
        if DAOIndex == 1
            stf(i).propVMAT.DAOAngleBorders = ([pln.propStf.DAOGantryAngles(DAOIndex) pln.propStf.DAOGantryAngles(DAOIndex+1)]+pln.propStf.DAOGantryAngles(DAOIndex))/2;
            
            lastDAOIndex = i;
            nextDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex+1));
            
            stf(i).propVMAT.lastDAOIndex = i;
            stf(i).propVMAT.nextDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex+1));
        elseif DAOIndex == length(pln.propStf.DAOGantryAngles)
            stf(i).propVMAT.DAOAngleBorders = ([pln.propStf.DAOGantryAngles(DAOIndex-1) pln.propStf.DAOGantryAngles(DAOIndex)]+pln.propStf.DAOGantryAngles(DAOIndex))/2;
            
            stf(i).propVMAT.lastDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex-1));
            stf(i).propVMAT.nextDAOIndex = i;
        else
            stf(i).propVMAT.DAOAngleBorders = ([pln.propStf.DAOGantryAngles(DAOIndex-1) pln.propStf.DAOGantryAngles(DAOIndex+1)]+pln.propStf.DAOGantryAngles(DAOIndex))/2;
            
            lastDAOIndex = i;
            nextDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex+1));
            
            stf(i).propVMAT.lastDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex-1));
            stf(i).propVMAT.nextDAOIndex = find(pln.propStf.gantryAngles == pln.propStf.DAOGantryAngles(DAOIndex+1));
        end
        stf(i).propVMAT.doseAngleDAO = ones(1,2);
        
        stf(i).propVMAT.DAOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.DAOAngleBorders(1) stf(i).propVMAT.DAOAngleBorders(2)-stf(i).gantryAngle];
        stf(i).propVMAT.DAOAngleBordersDiff = sum(stf(i).propVMAT.DAOAngleBorderCentreDiff);
        
        %This is the factor that relates the total time in the
        %optimized arc sector to the total time in the current dose
        %sector
        stf(i).propVMAT.timeFacCurr =  stf(i).propVMAT.doseAngleBordersDiff./stf(i).propVMAT.DAOAngleBordersDiff;
        
        %These are the factors that relate the total time in the
        %optimized arc sector to the total time in the previous and
        %next dose sectors
        stf(i).propVMAT.timeFac = zeros(1,2);
        
        if i == 1
            stf(i).propVMAT.timeFac(1) = 0;
            stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.DAOAngleBorderCentreDiff(2)/stf(i).propVMAT.DAOAngleBordersDiff;
        elseif i == length(pln.propStf.gantryAngles)
            stf(i).propVMAT.timeFac(1) = stf(i).propVMAT.DAOAngleBorderCentreDiff(1)/stf(i).propVMAT.DAOAngleBordersDiff;
            stf(i).propVMAT.timeFac(2) = 0;
        else
            stf(i).propVMAT.timeFac(1) = stf(i).propVMAT.DAOAngleBorderCentreDiff(1)/stf(i).propVMAT.DAOAngleBordersDiff;
            stf(i).propVMAT.timeFac(2) = stf(i).propVMAT.DAOAngleBorderCentreDiff(2)/stf(i).propVMAT.DAOAngleBordersDiff;
        end
        
    else
        if ~isfield(stf(stf(i).propVMAT.beamParentIndex).propVMAT,'beamSubChildrenGantryAngles') || isempty(stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles)
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = 0;
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles = nan(1000,1);
            stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex = nan(1000,1);
        end
        
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren = stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren+1;
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenGantryAngles(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = pln.propStf.gantryAngles(i);
        stf(stf(i).propVMAT.beamParentIndex).propVMAT.beamSubChildrenIndex(stf(stf(i).propVMAT.beamParentIndex).propVMAT.numOfBeamSubChildren) = i;
        
        stf(i).propVMAT.fracFromLastDAO = (pln.propStf.gantryAngles(nextDAOIndex)-pln.propStf.gantryAngles(i))./(pln.propStf.gantryAngles(nextDAOIndex)-pln.propStf.gantryAngles(lastDAOIndex));
        stf(i).propVMAT.lastDAOIndex = lastDAOIndex;
        stf(i).propVMAT.nextDAOIndex = nextDAOIndex;
    end
    
    
    if stf(i).propVMAT.FMOBeam
        % FMOAngleBorders are the angular borders over which an optimized
        % control point has influence
        FMOIndex = find(pln.propStf.FMOGantryAngles == pln.propStf.gantryAngles(i));
        
        if FMOIndex == 1
            
            stf(i).propVMAT.FMOAngleBorders = [min(pln.propStf.FMOGantryAngles(FMOIndex),pln.propStf.gantryAngles(1)) (pln.propStf.FMOGantryAngles(FMOIndex+1)+pln.propStf.FMOGantryAngles(FMOIndex))/2];
        elseif FMOIndex == length(pln.propStf.FMOGantryAngles)
            
            stf(i).propVMAT.FMOAngleBorders = [(pln.propStf.FMOGantryAngles(FMOIndex-1)+pln.propStf.FMOGantryAngles(FMOIndex))/2 max(pln.propStf.FMOGantryAngles(FMOIndex),pln.propStf.gantryAngles(end))];
        else
            
            stf(i).propVMAT.FMOAngleBorders = ([pln.propStf.FMOGantryAngles(FMOIndex-1) pln.propStf.FMOGantryAngles(FMOIndex+1)]+pln.propStf.FMOGantryAngles(FMOIndex))/2;
        end
        stf(i).propVMAT.FMOAngleBorderCentreDiff = [stf(i).gantryAngle-stf(i).propVMAT.FMOAngleBorders(1) stf(i).propVMAT.FMOAngleBorders(2)-stf(i).gantryAngle];
        stf(i).propVMAT.FMOAngleBordersDiff = sum(stf(i).propVMAT.FMOAngleBorderCentreDiff);
    end
    
    %% transformation of union of rays
    currMasterRayPosBEV = masterRayPosBEV;
    currMasterTargetPointBEV = masterTargetPointBEV;
    
    currMasterRayPosBEV(isnan(currMasterRayPosBEV)) = [];
    currMasterTargetPointBEV(isnan(currMasterTargetPointBEV)) = [];
    currMasterRayPosBEV = reshape(currMasterRayPosBEV,[],3);
    currMasterTargetPointBEV = reshape(currMasterTargetPointBEV,[],3);
    
    stf(i).numOfRays = size(currMasterRayPosBEV,1);
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    
    
    % source position in bev
    stf(i).sourcePoint_bev = [0 -SAD 0];
    
    % get (active) rotation matrix
    % transpose matrix because we are working with row vectors
    rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(i),pln.propStf.couchAngles(i)));
    
    
    stf(i).sourcePoint = stf(i).sourcePoint_bev*rotMat_vectors_T;
    
    % Save ray and target position in lps system.
    for j = 1:stf(i).numOfRays
        stf(i).ray(j).rayPos_bev = currMasterRayPosBEV(j,:);
        stf(i).ray(j).targetPoint_bev = currMasterTargetPointBEV(j,:);
        
        stf(i).ray(j).rayPos      = stf(i).ray(j).rayPos_bev*rotMat_vectors_T;
        stf(i).ray(j).targetPoint = stf(i).ray(j).targetPoint_bev*rotMat_vectors_T;
        if strcmp(pln.radiationMode,'photons')
            stf(i).ray(j).rayCorners_SCD = (repmat([0, machine.meta.SCD - SAD, 0],4,1)+ (machine.meta.SCD/SAD) * ...
                [currMasterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                currMasterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,+stf(i).bixelWidth/2];...
                currMasterRayPosBEV(j,:) + [-stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2];...
                currMasterRayPosBEV(j,:) + [+stf(i).bixelWidth/2,0,-stf(i).bixelWidth/2]])*rotMat_vectors_T;
        end
    end
    
    % loop over all rays to determine meta information for each ray
    stf(i).totalNumOfBixels = sum(stf(i).numOfBixelsPerRay);
    stf(i).numOfBixelsPerRay = ones(1,stf(i).numOfRays);
    
    for j = stf(i).numOfRays:-1:1
        
        % find appropriate energies for particles
        if strcmp(stf(i).radiationMode,'photons')
            
            % book keeping for photons
            stf(i).ray(j).energy = machine.data.energy;
        else
            error('Error generating stf struct: invalid radiation modality for VMAT.');
        end
    end
    
    matRad_progress(i,length(pln.propStf.gantryAngles));
end


%% final cleanup and calculation of factors we couldn't calc before
fprintf('matRad: VMAT post-processing (2/2)... ');

for i = 1:length(pln.propStf.gantryAngles)
    if stf(i).propVMAT.FMOBeam
        %remove NaNs from beamChildren and beamSubChildren
        if isfield(stf(i).propVMAT,'beamChildrenGantryAngles')
            stf(i).propVMAT.beamChildrenGantryAngles(isnan(stf(i).propVMAT.beamChildrenGantryAngles)) = [];
            stf(i).propVMAT.beamChildrenIndex(isnan(stf(i).propVMAT.beamChildrenIndex)) = [];
        else
            stf(i).propVMAT.numOfBeamChildren = 0;
        end
        if isfield(stf(i).propVMAT,'beamSubChildrenGantryAngles')
            stf(i).propVMAT.beamSubChildrenGantryAngles(isnan(stf(i).propVMAT.beamSubChildrenGantryAngles)) = [];
            stf(i).propVMAT.beamSubChildrenIndex(isnan(stf(i).propVMAT.beamSubChildrenIndex)) = [];
        else
            stf(i).propVMAT.numOfBeamSubChildren = 0;
        end
    end
    
    if ~stf(i).propVMAT.FMOBeam && ~stf(i).propVMAT.DAOBeam
        
        % for time interpolation
        stf(i).propVMAT.timeFracFromLastDAO = (stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2)-stf(i).propVMAT.doseAngleBorders(1))./stf(i).propVMAT.doseAngleBordersDiff;
        stf(i).propVMAT.timeFracFromNextDAO = (stf(i).propVMAT.doseAngleBorders(2)-stf(stf(i).propVMAT.lastDAOIndex).propVMAT.DAOAngleBorders(2))./stf(i).propVMAT.doseAngleBordersDiff;
        if stf(i).propVMAT.timeFracFromLastDAO > 1
            stf(i).propVMAT.timeFracFromLastDAO = 1;
        elseif stf(i).propVMAT.timeFracFromLastDAO < 0
            stf(i).propVMAT.timeFracFromLastDAO = 0;
        end
        if stf(i).propVMAT.timeFracFromNextDAO > 1
            stf(i).propVMAT.timeFracFromNextDAO = 1;
        elseif stf(i).propVMAT.timeFracFromNextDAO < 0
            stf(i).propVMAT.timeFracFromNextDAO = 0;
        end
    end
    
    matRad_progress(i,length(pln.propStf.gantryAngles));
end

