function updatedInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to translate the vector representation of the aperture
% shape and weight into an aperture info struct. At the same time, the
% updated bixel weight vector w is computed and a vector listing the
% correspondence between leaf tips and bixel indices for gradient
% calculation
%
% call
%   [updatedInfo,w,indVect] = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVect)
%
% input
%   apertureInfo:     aperture shape info struct
%   apertureInfoVect: aperture weights and shapes parameterized as vector
%
% output
%   updatedInfo: updated aperture shape info struct according to apertureInfoVect
%
% References
%
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

% function to update the apertureInfo struct after the each iteraton of the
% optimization

% initializing variables
updatedInfo = apertureInfo;

updatedInfo.apertureVector = apertureInfoVect;

if ~updatedInfo.runVMAT
    shapeInd = 1;
    
    indVect = NaN*ones(2*apertureInfo.doseTotalNumOfLeafPairs,1);
    offset = 0;
    
    % helper function to cope with numerical instabilities through rounding
    round2 = @(a,b) round(a*10^b)/10^b;
else
    
    % options for bixel and Jacobian calculation
    mlcOptions.bixelWidth = apertureInfo.bixelWidth;
    calcOptions.continuousAperture = updatedInfo.propVMAT.continuousAperture;
    vectorIndices.totalNumOfShapes = apertureInfo.totalNumOfShapes;
end

w = zeros(apertureInfo.totalNumOfBixels,1);

if updatedInfo.runVMAT && ~all([updatedInfo.propVMAT.beam.DAOBeam])
    j = 1;
    for i = 1:numel(updatedInfo.beam)
        if updatedInfo.propVMAT.beam(i).DAOBeam
            % update the shape weight
            % rescale the weight from the vector using the previous
            % iteration scaling factor
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset)./updatedInfo.beam(i).shape(j).jacobiScale;
            
            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
            updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
            updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU./updatedInfo.beam(i).time;
        end
    end
end

if updatedInfo.runVMAT
    %% ONLY SUPPORTED IN VMAT FOR NOW
    
    % Jacobian matrix to be used in the DAO gradient function
    % this tells us the gradient of a particular bixel with respect to an
    % element in the apertureVector (aperture weight or leaf position)
    % store as a vector for now, convert to sparse matrix later
    
    optBixelFactor = 7;
    % For optimized beams: 5 = (1 from weights) + (3 from left leaf positions (I, M, and F)) + (3 from
    % right leaf positions (I, M, and F))
    
    if updatedInfo.runVMAT
        intBixelFactor = 2*optBixelFactor+2;
        % For interpolated beams: multiply this number times 2 (influenced by the
        % one before and the one after), then add 2 (influenced by the time of the
        % times before and after)
    else
        intBixelFactor = 2*optBixelFactor;
        % For interpolated beams: multiply this number times 2 (influenced by the
        % one before and the one after)
    end
    
    % for the time (probability) gradients
    optBixelFactor = optBixelFactor+apertureInfo.totalNumOfShapes;
    intBixelFactor = intBixelFactor+apertureInfo.totalNumOfShapes;
    
    bixelJApVec_sz = (updatedInfo.totalNumOfOptBixels*optBixelFactor+(updatedInfo.totalNumOfBixels-updatedInfo.totalNumOfOptBixels)*intBixelFactor)*2;
    
    bixelJApVec_vec = zeros(1,bixelJApVec_sz);
    
    % vector indices
    bixelJApVec_i = nan(1,bixelJApVec_sz);
    % bixel indices
    bixelJApVec_j = zeros(1,bixelJApVec_sz);
    % offset
    bixelJApVec_offset = 0;
end

%% update the shapeMaps
% here the new colimator positions are used to create new shapeMaps that
% now include decimal values instead of binary

calcOptions.saveJacobian = true;

% loop over all beams
for i = 1:numel(updatedInfo.beam)
    
    %posOfRightCornerPixel = apertureInfo.beam(i).posOfCornerBixel(1) + (size(apertureInfo.beam(i).bixelIndMap,2)-1)*apertureInfo.bixelWidth;
    
    % pre compute left and right bixel edges
    edges_l = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ((1:size(apertureInfo.beam(i).bixelIndMap,2))-1-1/2)*updatedInfo.bixelWidth;
    edges_r = updatedInfo.beam(i).posOfCornerBixel(1)...
        + ((1:size(apertureInfo.beam(i).bixelIndMap,2))-1+1/2)*updatedInfo.bixelWidth;
    
    % get dimensions of 2d matrices that store shape/bixel information
    n = apertureInfo.beam(i).numOfActiveLeafPairs;
    
    % loop over all shapes
    if updatedInfo.runVMAT
        numOfShapes = 1;
        calcOptions.DAOBeam = updatedInfo.propVMAT.beam(i).DAOBeam;
    else
        numOfShapes = updatedInfo.beam(i).numOfShapes;
    end
    
    mlcOptions.lim_l = apertureInfo.beam(i).lim_l;
    mlcOptions.lim_r = apertureInfo.beam(i).lim_r;
    mlcOptions.edges_l = edges_l;
    mlcOptions.edges_r = edges_r;
    mlcOptions.centres = (edges_l+edges_r)/2;
    mlcOptions.widths = edges_r-edges_l;
    mlcOptions.n = n;
    mlcOptions.numBix = size(apertureInfo.beam(i).bixelIndMap,2);
    mlcOptions.bixelIndMap = apertureInfo.beam(i).bixelIndMap;
    
    for j = 1:numOfShapes
        
        if ~updatedInfo.runVMAT || updatedInfo.propVMAT.beam(i).DAOBeam
            % either this is not VMAT, or if it is VMAT, this is a DAO beam
            
            % update the shape weight
            updatedInfo.beam(i).shape(j).weight = apertureInfoVect(updatedInfo.beam(i).shape(j).weightOffset)./updatedInfo.beam(i).shape(j).jacobiScale;
            
            if updatedInfo.runVMAT
                updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).weight*updatedInfo.weightToMU;
                updatedInfo.beam(i).time = apertureInfoVect((updatedInfo.totalNumOfShapes+updatedInfo.totalNumOfLeafPairs*2)+updatedInfo.propVMAT.beam(i).DAOIndex)*updatedInfo.propVMAT.beam(i).timeFacCurr;
                updatedInfo.beam(i).gantryRot = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff/updatedInfo.beam(i).time;
                updatedInfo.beam(i).shape(j).MURate = updatedInfo.beam(i).shape(j).MU./updatedInfo.beam(i).time;
            end
            
            if ~updatedInfo.runVMAT || ~updatedInfo.propVMAT.continuousAperture
                % extract left and right leaf positions from shape vector
                vectorIx_L = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n)-1);
                vectorIx_R = vectorIx_L+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos  = apertureInfoVect(vectorIx_L);
                rightLeafPos = apertureInfoVect(vectorIx_R);
                
                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_I = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_F = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos;
            else
                % extract left and right leaf positions from shape vector
                vectorIx_LI = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n)-1);
                vectorIx_RI = vectorIx_LI+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_I = apertureInfoVect(vectorIx_LI);
                rightLeafPos_I = apertureInfoVect(vectorIx_RI);
                
                vectorIx_LF = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n)-1);
                vectorIx_RF = vectorIx_LF+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_F = apertureInfoVect(vectorIx_LF);
                rightLeafPos_F = apertureInfoVect(vectorIx_RF);
                
                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos_I  = leftLeafPos_I;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos_I;
                
                updatedInfo.beam(i).shape(j).leftLeafPos_F  = leftLeafPos_F;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos_F;
            end
            
        else
            % this is an interpolated beam
            
            %MURate is interpolated between MURates of optimized apertures
            updatedInfo.beam(i).gantryRot = 1./(updatedInfo.propVMAT.beam(i).timeFracFromLastDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).gantryRot+updatedInfo.propVMAT.beam(i).timeFracFromNextDAO./updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).gantryRot);
            updatedInfo.beam(i).time = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff./updatedInfo.beam(i).gantryRot;
            updatedInfo.beam(i).shape(j).MURate = updatedInfo.propVMAT.beam(i).fracFromLastDAO*updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).MURate+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).MURate;
            
            % calculate MU, weight
            updatedInfo.beam(i).shape(j).MU = updatedInfo.beam(i).shape(j).MURate.*updatedInfo.beam(i).time;
            updatedInfo.beam(i).shape(j).weight = updatedInfo.beam(i).shape(j).MU./updatedInfo.weightToMU;
            
            if ~updatedInfo.propVMAT.continuousAperture
                
                fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO*ones(n,1);
                fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO*ones(n,1);
                fracFromNextOptI = (1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*ones(n,1);
                fracFromNextOptF = (1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*ones(n,1);
                
                % obtain leaf positions at last DAO beam
                vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset + ((1:n)-1);
                vectorIx_RF_last = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_last = apertureInfoVect(vectorIx_LF_last);
                rightLeafPos_last = apertureInfoVect(vectorIx_RF_last);
                
                % obtain leaf positions at next DAO beam
                vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset + ((1:n)-1);
                vectorIx_RI_next = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_next = apertureInfoVect(vectorIx_LI_next);
                rightLeafPos_next = apertureInfoVect(vectorIx_RI_next);
                
                % interpolate leaf positions
                leftLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*leftLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*leftLeafPos_next;
                rightLeafPos = updatedInfo.propVMAT.beam(i).fracFromLastDAO*rightLeafPos_last+(1-updatedInfo.propVMAT.beam(i).fracFromLastDAO)*rightLeafPos_next;
                
                % update information in shape structure
                updatedInfo.beam(i).shape(j).leftLeafPos  = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_I = leftLeafPos;
                updatedInfo.beam(i).shape(j).leftLeafPos_F = leftLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = rightLeafPos;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = rightLeafPos;
            else
                
                fracFromLastOpt = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                fracFromLastOptI = updatedInfo.propVMAT.beam(i).fracFromLastDAO_I*ones(n,1);
                fracFromLastOptF = updatedInfo.propVMAT.beam(i).fracFromLastDAO_F*ones(n,1);
                fracFromNextOptI = updatedInfo.propVMAT.beam(i).fracFromNextDAO_I*ones(n,1);
                fracFromNextOptF = updatedInfo.propVMAT.beam(i).fracFromNextDAO_F*ones(n,1);
                
                % obtain leaf positions at last DAO beam
                vectorIx_LF_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset(2) + ((1:n)-1);
                vectorIx_RF_last = vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_F_last = apertureInfoVect(vectorIx_LF_last);
                rightLeafPos_F_last = apertureInfoVect(vectorIx_RF_last);
                
                % obtain leaf positions at next DAO beam
                vectorIx_LI_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset(1) + ((1:n)-1);
                vectorIx_RI_next = vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs;
                leftLeafPos_I_next = apertureInfoVect(vectorIx_LI_next);
                rightLeafPos_I_next = apertureInfoVect(vectorIx_RI_next);
                
                % interpolate leaf positions
                updatedInfo.beam(i).shape(j).leftLeafPos_I = fracFromLastOptI.*leftLeafPos_F_last+fracFromNextOptI.*leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPos_I = fracFromLastOptI.*rightLeafPos_F_last+fracFromNextOptI.*rightLeafPos_I_next;
                
                updatedInfo.beam(i).shape(j).leftLeafPos_F = fracFromLastOptF.*leftLeafPos_F_last+fracFromNextOptF.*leftLeafPos_I_next;
                updatedInfo.beam(i).shape(j).rightLeafPos_F = fracFromLastOptF.*rightLeafPos_F_last+fracFromNextOptF.*rightLeafPos_I_next;
            end
        end
        
        if ~updatedInfo.runVMAT
            
            % rounding for numerical stability
            leftLeafPos  = round2(leftLeafPos,10);
            rightLeafPos = round2(rightLeafPos,10);
            
            % check overshoot of leaf positions
            leftLeafPos(leftLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(leftLeafPos <= apertureInfo.beam(i).lim_l);
            rightLeafPos(rightLeafPos <= apertureInfo.beam(i).lim_l) = apertureInfo.beam(i).lim_l(rightLeafPos <= apertureInfo.beam(i).lim_l);
            leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(leftLeafPos >= apertureInfo.beam(i).lim_r);
            rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) = apertureInfo.beam(i).lim_r(rightLeafPos >= apertureInfo.beam(i).lim_r);
            
            %
            xPosIndLeftLeaf  = round((leftLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            xPosIndRightLeaf = round((rightLeafPos - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            
            %
            xPosIndLeftLeaf_lim  = floor((apertureInfo.beam(i).lim_l - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth+1);
            xPosIndRightLeaf_lim = ceil((apertureInfo.beam(i).lim_r - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth + 1);
            
            xPosIndLeftLeaf(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim) = xPosIndLeftLeaf_lim(xPosIndLeftLeaf <= xPosIndLeftLeaf_lim)+1;
            xPosIndRightLeaf(xPosIndRightLeaf >= xPosIndRightLeaf_lim) = xPosIndRightLeaf_lim(xPosIndRightLeaf >= xPosIndRightLeaf_lim)-1;
            
            % check limits because of rounding off issues at maximum, i.e.,
            % enforce round(X.5) -> X
            % LeafPos can occasionally go slightly beyond lim_r, so changed
            % == check to >=
            xPosIndLeftLeaf(leftLeafPos >= apertureInfo.beam(i).lim_r) = round(...
                .5 + (leftLeafPos(leftLeafPos >= apertureInfo.beam(i).lim_r) ...
                - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
            
            xPosIndRightLeaf(rightLeafPos >= apertureInfo.beam(i).lim_r) = round(...
                .5 + (rightLeafPos(rightLeafPos >= apertureInfo.beam(i).lim_r) ...
                - apertureInfo.beam(i).posOfCornerBixel(1))/apertureInfo.bixelWidth);
            
            % find the bixel index that the leaves currently touch
            bixelIndLeftLeaf  = apertureInfo.beam(i).bixelIndMap((xPosIndLeftLeaf-1)*n+[1:n]');
            bixelIndRightLeaf = apertureInfo.beam(i).bixelIndMap((xPosIndRightLeaf-1)*n+[1:n]');
            
            if any(isnan(bixelIndLeftLeaf)) || any(isnan(bixelIndRightLeaf))
                error('cannot map leaf position to bixel index');
            end
            
            % store information in index vector for gradient calculation
            indVect(offset+(1:n)) = bixelIndLeftLeaf;
            indVect(offset+(1:n)+apertureInfo.doseTotalNumOfLeafPairs) = bixelIndRightLeaf;
            offset = offset+n;
            
            % calculate opening fraction for every bixel in shape to construct
            % bixel weight vector
            
            coveredByLeftLeaf  = bsxfun(@minus,leftLeafPos,edges_l)  / updatedInfo.bixelWidth;
            coveredByRightLeaf = bsxfun(@minus,edges_r,rightLeafPos) / updatedInfo.bixelWidth;
            
            tempMap = 1 - (coveredByLeftLeaf  + abs(coveredByLeftLeaf))  / 2 ...
                - (coveredByRightLeaf + abs(coveredByRightLeaf)) / 2;
            
            % find open bixels
            tempMapIx = tempMap > 0;
            
            currBixelIx = apertureInfo.beam(i).bixelIndMap(tempMapIx);
            w(currBixelIx) = w(currBixelIx) + tempMap(tempMapIx)*updatedInfo.beam(i).shape(j).weight;
            
            % save the tempMap (we need to apply a positivity operator !)
            updatedInfo.beam(i).shape(j).shapeMap = (tempMap  + abs(tempMap))  / 2;
            
            % increment shape index
            shapeInd = shapeInd +1;
        end
        
    end
    
    if updatedInfo.runVMAT
        
        for j = 1:numOfShapes
            
            % shapeMap
            shapeMap = zeros(size(updatedInfo.beam(i).bixelIndMap));
            % sumGradSq
            sumGradSq = 0;
            
            % insert variables
            vectorIndices.tIx_Vec       = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+(1:apertureInfo.totalNumOfShapes);
            
            variables.weight            = updatedInfo.beam(i).shape(j).weight;
            variables.leftLeafPos_I     = updatedInfo.beam(i).shape(j).leftLeafPos_I;
            variables.leftLeafPos_F     = updatedInfo.beam(i).shape(j).leftLeafPos_F;
            variables.rightLeafPos_I    = updatedInfo.beam(i).shape(j).rightLeafPos_I;
            variables.rightLeafPos_F    = updatedInfo.beam(i).shape(j).rightLeafPos_F;
            
            if updatedInfo.propVMAT.beam(i).DAOBeam
                
                variables.jacobiScale = updatedInfo.beam(i).shape(1).jacobiScale;
                
                vectorIndices.DAOindex      = updatedInfo.propVMAT.beam(i).DAOIndex;
                if updatedInfo.propVMAT.continuousAperture
                    vectorIndices.vectorIx_LI   = updatedInfo.beam(i).shape(j).vectorOffset(1) + ((1:n)-1);
                    vectorIndices.vectorIx_LF   = updatedInfo.beam(i).shape(j).vectorOffset(2) + ((1:n)-1);
                    vectorIndices.vectorIx_RI   = vectorIndices.vectorIx_LI+apertureInfo.totalNumOfLeafPairs;
                    vectorIndices.vectorIx_RF   = vectorIndices.vectorIx_LF+apertureInfo.totalNumOfLeafPairs;
                else
                    vectorIndices.vectorIx_LI   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n)-1);
                    vectorIndices.vectorIx_LF   = updatedInfo.beam(i).shape(j).vectorOffset + ((1:n)-1);
                    vectorIndices.vectorIx_RI   = vectorIndices.vectorIx_LI+apertureInfo.totalNumOfLeafPairs;
                    vectorIndices.vectorIx_RF   = vectorIndices.vectorIx_LF+apertureInfo.totalNumOfLeafPairs;
                end
            else
                
                variables.weight_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).weight;
                variables.weight_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).weight;
                
                variables.jacobiScale_last    = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(1).jacobiScale;
                variables.jacobiScale_next    = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(1).jacobiScale;
                
                variables.time_last = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).time;
                variables.time_next = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).time;
                variables.time      = updatedInfo.beam(i).time;
                
                variables.fracFromLastOptI  = fracFromLastOptI;
                variables.fracFromLastOptF  = fracFromLastOptF;
                variables.fracFromNextOptI  = fracFromNextOptI;
                variables.fracFromNextOptF  = fracFromNextOptF;
                variables.fracFromLastOpt   = fracFromLastOpt;
                
                variables.doseAngleBordersDiff      = updatedInfo.propVMAT.beam(i).doseAngleBordersDiff;
                variables.doseAngleBordersDiff_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff;
                variables.doseAngleBordersDiff_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff;
                variables.timeFacCurr_last          = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr;
                variables.timeFacCurr_next          = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr;
                variables.fracFromLastDAO           = updatedInfo.propVMAT.beam(i).fracFromLastDAO;
                variables.timeFracFromLastDAO       = updatedInfo.propVMAT.beam(i).timeFracFromLastDAO;
                variables.timeFracFromNextDAO       = updatedInfo.propVMAT.beam(i).timeFracFromNextDAO;
                
                vectorIndices.DAOindex_last = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
                vectorIndices.DAOindex_next = updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
                vectorIndices.tIx_last      = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).DAOIndex;
                vectorIndices.tIx_next      = (apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+updatedInfo.propVMAT.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).DAOIndex;
                
                if updatedInfo.propVMAT.continuousAperture
                    vectorIndices.vectorIx_LF_last  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset(2) + ((1:n)-1);
                    vectorIndices.vectorIx_LI_next  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset(1) + ((1:n)-1);
                    vectorIndices.vectorIx_RF_last  = vectorIndices.vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs;
                    vectorIndices.vectorIx_RI_next  = vectorIndices.vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs;
                else
                    vectorIndices.vectorIx_LF_last  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).lastDAOIndex).shape(j).vectorOffset + ((1:n)-1);
                    vectorIndices.vectorIx_LI_next  = updatedInfo.beam(updatedInfo.propVMAT.beam(i).nextDAOIndex).shape(j).vectorOffset + ((1:n)-1);
                    vectorIndices.vectorIx_RF_last  = vectorIndices.vectorIx_LF_last+apertureInfo.totalNumOfLeafPairs;
                    vectorIndices.vectorIx_RI_next  = vectorIndices.vectorIx_LI_next+apertureInfo.totalNumOfLeafPairs;
                end
            end
            
            counters.bixelJApVec_offset = bixelJApVec_offset;
            
            % calculate bixel weight and derivative in function
            [w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap,counters] = ...
                matRad_bixWeightAndGrad(calcOptions,mlcOptions,variables,vectorIndices,counters,w,bixelJApVec_vec,bixelJApVec_i,bixelJApVec_j,sumGradSq,shapeMap);
            
            bixelJApVec_offset = counters.bixelJApVec_offset;
            
            % update shapeMap
            updatedInfo.beam(i).shape(j).shapeMap = shapeMap;
            % update sumGradSq
            % FIX THIS FOR INTERPOLATED ANGLES???
            updatedInfo.beam(i).shape(j).sumGradSq = sumGradSq;
            
        end
    end
end


% save bixelWeight, apertureVector
updatedInfo.bixelWeights = w;
updatedInfo.apertureVector = apertureInfoVect;

if updatedInfo.runVMAT
    % save Jacobian between bixelWeight, apertureVector
    
    deleteInd_i = isnan(bixelJApVec_i);
    deleteInd_j = bixelJApVec_j == 0;
    if ~all(deleteInd_i == deleteInd_j)
        error('Jacobian deletion mismatch');
    else
        bixelJApVec_i(deleteInd_i) = [];
        bixelJApVec_j(deleteInd_i) = [];
        bixelJApVec_vec(deleteInd_i) = [];
    end
    updatedInfo.bixelJApVec = sparse(bixelJApVec_i,bixelJApVec_j,bixelJApVec_vec,numel(apertureInfoVect),updatedInfo.totalNumOfBixels);
else
    % save indVect
    updatedInfo.bixelIndices = indVect;
end

end