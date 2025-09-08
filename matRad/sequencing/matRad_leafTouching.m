function apertureInfo = matRad_leafTouching(apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to improve instances of leaf touching by moving leaves
% from the centre to sweep with the non-touching leaves.
% 
% Currently only works with VMAT, add option to work with IMRT (not as
% crucial)
%
% call
%   apertureInfo = matRad_leafTouching(apertureInfo)
%
% input
%   apertureInfo: matRad aperture weight and shape info struct
%
% output
%   apertureInfo: matRad aperture weight and shape info struct
%
% References
%   
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


%initialize
dimZ = apertureInfo.beam(1).numOfActiveLeafPairs;
numBeams = nnz([apertureInfo.propVMAT.beam.DAOBeam]);
if ~isfield(apertureInfo.beam(1).shape(1),'leftLeafPos_I')
    % Each non-interpolated beam should have 1 left/right leaf position
    leftLeafPoss = nan(dimZ,numBeams);
    rightLeafPoss = nan(dimZ,numBeams);
    gantryAngles = zeros(1,numBeams);
else
    % Each non-interpolated beam should have 2 left/right leaf positions
    leftLeafPoss = nan(dimZ,2*numBeams);
    rightLeafPoss = nan(dimZ,2*numBeams);
    gantryAngles = zeros(1,2*numBeams);
end
initBorderGantryAngles = unique([apertureInfo.propVMAT.beam.FMOAngleBorders]);
initBorderLeftLeafPoss = nan(dimZ,numel(initBorderGantryAngles));

l = 1;
m = 1;
%collect all leaf positions
for k = 1:numel(apertureInfo.beam)
    if (k ~= 1 && apertureInfo.beam(k).gantryAngle == apertureInfo.beam(k-1).gantryAngle) || ~apertureInfo.propVMAT.beam(k).DAOBeam
        continue
    end
    
    if ~isfield(apertureInfo.beam(1).shape(1),'leftLeafPos_I')
        leftLeafPoss(:,l) = apertureInfo.beam(k).shape(1).leftLeafPos;
        rightLeafPoss(:,l) = apertureInfo.beam(k).shape(1).rightLeafPos;
        gantryAngles(l) = apertureInfo.beam(k).gantryAngle;
        
        l = l+1;
    else
        leftLeafPoss(:,l) = apertureInfo.beam(k).shape(1).leftLeafPos_I;
        rightLeafPoss(:,l) = apertureInfo.beam(k).shape(1).rightLeafPos_I;
        gantryAngles(l) = apertureInfo.beam(k).doseAngleBorders(1);
        
        l = l+1;
        
        leftLeafPoss(:,l) = apertureInfo.beam(k).shape(1).leftLeafPos_F;
        rightLeafPoss(:,l) = apertureInfo.beam(k).shape(1).rightLeafPos_F;
        gantryAngles(l) = apertureInfo.beam(k).doseAngleBorders(2);
        
        l = l+1;
    end
    
    %Only important when cleaning up instances of opposing
    %leaves touching.
    if apertureInfo.propVMAT.beam(k).FMOBeam
        if apertureInfo.propVMAT.beam(k).leafDir == 1
            %This means that the current arc sector is moving
            %in the normal direction (L-R).
            initBorderLeftLeafPoss(:,m) = apertureInfo.beam(k).lim_l;
            
        elseif apertureInfo.propVMAT.beam(k).leafDir == -1
            %This means that the current arc sector is moving
            %in the reverse direction (R-L).
            initBorderLeftLeafPoss(:,m) = apertureInfo.beam(k).lim_r;
        end
        m = m+1;
        
        %end of last sector
        if m == numel(initBorderGantryAngles)
            %This gives ending angle of the current sector.
            if apertureInfo.propVMAT.beam(k).leafDir == 1
                %This means that the current arc sector is moving
                %in the normal direction (L-R), so the next arc
                %sector is moving opposite
                initBorderLeftLeafPoss(:,m) = apertureInfo.beam(k).lim_r;
            elseif apertureInfo.propVMAT.beam(k).leafDir == -1
                %This means that the current arc sector is moving
                %in the reverse direction (R-L), so the next
                %arc sector is moving opposite
                initBorderLeftLeafPoss(:,m) = apertureInfo.beam(k).lim_l;
            end
        end
    end
end

[gantryAngles,ind] = unique(gantryAngles);
leftLeafPoss = leftLeafPoss(:,ind);
rightLeafPoss = rightLeafPoss(:,ind);

%Any time leaf pairs are touching, they are set to
%be in the middle of the field.  Instead, move them
%so that they are still touching, but that they
%follow the motion of the MLCs across the field.
for row = 1:dimZ
    
    touchingInd = find(leftLeafPoss(row,:) == rightLeafPoss(row,:));
    
    if ~exist('leftLeafPossAug','var')
        %leftLeafPossAug = [reshape(mean([leftLeafPoss(:) rightLeafPoss(:)],2),size(leftLeafPoss)),borderLeftLeafPoss];
        leftLeafPossAugTemp = reshape(mean([leftLeafPoss(:) rightLeafPoss(:)],2),size(leftLeafPoss));
        
        numRep = 0;
        repInd = nan(size(gantryAngles));
        for j = 1:numel(gantryAngles)
            if any(gantryAngles(j) == initBorderGantryAngles)
                %replace leaf positions with the ones at
                %the borders (eliminates repetitions)
                numRep = numRep+1;
                %these are the gantry angles that are
                %repeated
                repInd(numRep) = j;
                
                delInd = find(gantryAngles(j) == initBorderGantryAngles);
                leftLeafPossAugTemp(:,j) = initBorderLeftLeafPoss(:,delInd);
                initBorderLeftLeafPoss(:,delInd) = [];
                initBorderGantryAngles(delInd) = [];
            end
        end
        repInd(isnan(repInd)) = [];
        leftLeafPossAug = [leftLeafPossAugTemp,initBorderLeftLeafPoss];
        gantryAnglesAug = [gantryAngles,initBorderGantryAngles];
    end
    notTouchingInd = [setdiff(1:numBeams,touchingInd),repInd];
    notTouchingInd = unique(notTouchingInd);
    %make sure to include the repeated ones in the
    %interpolation!
    
    notTouchingIndAug = [notTouchingInd,(1+numel(gantryAngles)):(numel(gantryAngles)+numel(initBorderGantryAngles))];
    
    leftLeafPoss(row,touchingInd) = interp1(gantryAnglesAug(notTouchingIndAug),leftLeafPossAug(row,notTouchingIndAug),gantryAngles(touchingInd))-0.5;
    rightLeafPoss(row,touchingInd) = leftLeafPoss(row,touchingInd)+1;
end


%finally, set new leaf positions
for i = 1:numel(apertureInfo.beam)
    apertureInfo.beam(i).shape(1).leftLeafPos = max((interp1(gantryAngles',leftLeafPoss',apertureInfo.beam(i).gantryAngle))',apertureInfo.beam(i).lim_l);
    apertureInfo.beam(i).shape(1).rightLeafPos = min((interp1(gantryAngles',rightLeafPoss',apertureInfo.beam(i).gantryAngle))',apertureInfo.beam(i).lim_r);
    
    apertureInfo.beam(i).shape(1).leftLeafPos_I = max((interp1(gantryAngles',leftLeafPoss',apertureInfo.propVMAT.beam(i).doseAngleBorders(1)))',apertureInfo.beam(i).lim_l);
    apertureInfo.beam(i).shape(1).rightLeafPos_I = min((interp1(gantryAngles',rightLeafPoss',apertureInfo.propVMAT.beam(i).doseAngleBorders(1)))',apertureInfo.beam(i).lim_r);
    
    apertureInfo.beam(i).shape(1).leftLeafPos_F = max((interp1(gantryAngles',leftLeafPoss',apertureInfo.propVMAT.beam(i).doseAngleBorders(2)))',apertureInfo.beam(i).lim_l);
    apertureInfo.beam(i).shape(1).rightLeafPos_F = min((interp1(gantryAngles',rightLeafPoss',apertureInfo.propVMAT.beam(i).doseAngleBorders(2)))',apertureInfo.beam(i).lim_r);
end


end

