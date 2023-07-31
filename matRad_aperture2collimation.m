function [pln,stf] = matRad_aperture2collimation(pln,stf,sequencing,apertureInfo)
% matRad function to convert sequencing information / aperture information
% into collimation information in pln and stf for field-based dose
% calculation
% 
% call
%   [pln,stf] = matRad_aperture2collimation(pln,stf,sequencing,apertureInfo)
%
% input
%   pln:            pln file used to generate the sequenced plan
%   stf:            stf file used to generate the sequenced plan
%   sequencing:     sequencing information (from resultGUI)
%   apertureInfo:   apertureInfo (from resultGUI)
%
% output
%   pln:            matRad pln struct with collimation information
%   stf:            matRad stf struct with shapes instead of beamlets 
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% Author: wahln
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Debug visualization
visBool = false;

bixelWidth = apertureInfo.bixelWidth;
leafWidth = bixelWidth;
convResolution = 0.5; %[mm]

%The collimator limits are infered here from the apertureInfo. This could
%be handled differently by explicitly storing collimator info in the base
%data?
symmetricMLClimits = vertcat(apertureInfo.beam.MLCWindow);
symmetricMLClimits = max(abs(symmetricMLClimits));
fieldWidth = 2*max(symmetricMLClimits);

%modify basic pln variables
pln.propStf.bixelWidth = 'field';
pln.propStf.collimation.convResolution = 0.5; %[mm]
pln.propStf.collimation.fieldWidth = fieldWidth;
pln.propStf.collimation.leafWidth = leafWidth;

%
%[bixelFieldX,bixelFieldY] = ndgrid(-fieldWidth/2:bixelWidth:fieldWidth/2,-fieldWidth/2:leafWidth:fieldWidth/2);
[convFieldX,convFieldY] = meshgrid(-fieldWidth/2:convResolution:fieldWidth/2);

%TODO: Not used in calcPhotonDose but imported from DICOM
%pln.propStf.collimation.Devices ...
%pln.propStf.collimation.numOfFields
%pln.propStf.collimation.beamMeterset 

for iBeam = 1:numel(stf)
    stfTmp = stf(iBeam);
    beamSequencing = sequencing.beam(iBeam);
    beamAperture = apertureInfo.beam(iBeam);
    
    stfTmp.bixelWidth = 'field';
    
    nShapes = beamSequencing.numOfShapes;

    stfTmp.numOfRays = 1;%
    stfTmp.numOfBixelsPerRay = nShapes;
    stfTmp.totalNumOfBixels = nShapes;
    
    ray = struct();
    ray.rayPos_bev = [0 0 0];
    ray.targetPoint_bev = [0 stfTmp.SAD 0];
    ray.weight = 1;
    ray.energy = stfTmp.ray(1).energy;
    ray.beamletCornersAtIso = stfTmp.ray(1).beamletCornersAtIso;
    ray.rayCorners_SCD = stfTmp.ray(1).rayCorners_SCD;

    %ray.shape = beamSequencing.sum;
    shapeTotalF = zeros(size(convFieldX));

    ray.shapes = struct();
    for iShape = 1:nShapes
        currShape = beamAperture.shape(iShape);
        activeLeafPairPosY = beamAperture.leafPairPos;
        F = zeros(size(convFieldX));
        if visBool
            hF = figure; imagesc(F); title(sprintf('Beam %d, Shape %d',iBeam,iShape)); hold on;
        end
        for iLeafPair = 1:numel(activeLeafPairPosY)
            posY = activeLeafPairPosY(iLeafPair);
            ixY = convFieldY >= posY-leafWidth/2 & convFieldY < posY + leafWidth/2;
            ixX = convFieldX >= currShape.leftLeafPos(iLeafPair) & convFieldX < currShape.rightLeafPos(iLeafPair);
            ix = ixX & ixY;            
            F(ix) = 1;
            if visBool
                figure(hF); imagesc(F); drawnow; pause(0.1);
            end
        end

        if visBool
            pause(1); close(hF);
        end

        F = F*currShape.weight;
        shapeTotalF = shapeTotalF + F;

        ray.shapes(iShape).convFluence = F;
        ray.shapes(iShape).shapeMap = currShape.shapeMap;
        ray.shapes(iShape).weight = currShape.weight;
        ray.shapes(iShape).leftLeafPos = currShape.leftLeafPos;
        ray.shapes(iShape).rightLeafPos = currShape.rightLeafPos;
        ray.shapes(iShape).leafPairCenterPos = activeLeafPairPosY;
    end

    ray.shape = shapeTotalF;
    ray.weight = ones(1,nShapes);
    stfTmp.ray = ray;

    stf(iBeam) = stfTmp;
end


