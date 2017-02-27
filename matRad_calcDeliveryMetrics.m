function result = matRad_calcDeliveryMetrics(result,pln,stf)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad delivery metric calculation
% 
% call
%   matRad_calcDeliveryMetrics(result,pln)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   pln:                matRad plan meta information struct
%
% output
%   All plans: total MU
%   VMAT plans: total time, leaf speed, MU rate, and gantry rotation speed
%   distributions
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

apertureInfo = result.apertureInfo;

apertureInfo.planMU = 0;
apertureInfo.planArea = 0;
apertureInfo.planModulation = 0;

apertureMU = nan(1000,1);
apertureArea = nan(1000,1);
l = 1;
if pln.VMAT
    
    initInd = find([stf(:).initializeBeam]);
    
    for i = initInd
        apertureInfo.beam(i).beamMU = 0;
        apertureInfo.beam(i).beamArea = 0;
        apertureInfo.beam(i).unionShapeMap = 0*apertureInfo.beam(i).shape(1).shapeMap;
        
        for j = stf(i).beamChildrenIndex'
            apertureInfo.beam(j).shape(1).MU = apertureInfo.weightToMU*apertureInfo.beam(j).shape(1).weight;
            apertureInfo.beam(j).shape(1).apertureArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(j).shape(1).shapeMap(:));
            
            apertureInfo.beam(i).beamMU = apertureInfo.beam(i).beamMU+apertureInfo.beam(j).shape(1).MU;
            apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea+apertureInfo.beam(j).shape(1).MU*apertureInfo.beam(j).shape(1).apertureArea;
            
            apertureInfo.beam(i).unionShapeMap = max(apertureInfo.beam(i).unionShapeMap,apertureInfo.beam(j).shape(1).shapeMap);
            
            apertureMU(l) = apertureInfo.beam(j).shape(1).MU;
            apertureArea(l) = apertureInfo.beam(j).shape(1).apertureArea;
            l = l+1;
        end
        apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea./apertureInfo.beam(i).beamMU;
        apertureInfo.beam(i).unionArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(i).unionShapeMap(:));
        apertureInfo.beam(i).beamModulation = 1-apertureInfo.beam(i).beamArea./apertureInfo.beam(i).unionArea;
        apertureInfo.beam(i).beamK = apertureInfo.beam(i).numOfShapes*(1-apertureInfo.beam(i).beamModulation);
        
        apertureInfo.planMU = apertureInfo.planMU+apertureInfo.beam(i).beamMU;
        apertureInfo.planArea = apertureInfo.planArea+apertureInfo.beam(i).beamArea*apertureInfo.beam(i).beamMU;
        apertureInfo.planModulation = apertureInfo.planModulation+apertureInfo.beam(i).beamModulation*apertureInfo.beam(i).beamMU;
    end
else
    for i = 1:size(apertureInfo.beam,2)
        if apertureInfo.beam(i).numOfShapes ~= 0
            
            apertureInfo.beam(i).beamMU = 0;
            apertureInfo.beam(i).beamArea = 0;
            apertureInfo.beam(i).unionShapeMap = 0*apertureInfo.beam(i).shape(1).shapeMap;
            
            for j = 1:apertureInfo.beam(i).numOfShapes
                apertureInfo.beam(i).shape(j).MU = apertureInfo.weightToMU*apertureInfo.beam(i).shape(j).weight;
                apertureInfo.beam(i).shape(j).apertureArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(i).shape(j).shapeMap(:));
                
                apertureInfo.beam(i).beamMU = apertureInfo.beam(i).beamMU+apertureInfo.beam(i).shape(j).MU;
                apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea+apertureInfo.beam(i).shape(j).MU*apertureInfo.beam(i).shape(j).apertureArea;
                
                apertureInfo.beam(i).unionShapeMap = max(apertureInfo.beam(i).unionShapeMap,apertureInfo.beam(i).shape(j).shapeMap);
                
                apertureMU(l) = apertureInfo.beam(i).shape(j).MU;
                apertureArea(l) = apertureInfo.beam(i).shape(j).apertureArea;
                l = l+1;
            end
            apertureInfo.beam(i).beamArea = apertureInfo.beam(i).beamArea./apertureInfo.beam(i).beamMU;
            apertureInfo.beam(i).unionArea = (apertureInfo.bixelWidth/10)^2*sum(apertureInfo.beam(i).unionShapeMap(:));
            apertureInfo.beam(i).beamModulation = 1-apertureInfo.beam(i).beamArea./apertureInfo.beam(i).unionArea;
            apertureInfo.beam(i).beamK = apertureInfo.beam(i).numOfShapes*(1-apertureInfo.beam(i).beamModulation);
            
            apertureInfo.planMU = apertureInfo.planMU+apertureInfo.beam(i).beamMU;
            apertureInfo.planArea = apertureInfo.planArea+apertureInfo.beam(i).beamArea*apertureInfo.beam(i).beamMU;
            apertureInfo.planModulation = apertureInfo.planModulation+apertureInfo.beam(i).beamModulation*apertureInfo.beam(i).beamMU;
            
        end
    end
end

apertureInfo.planArea = apertureInfo.planArea./apertureInfo.planMU;
apertureInfo.planModulation = apertureInfo.planModulation./apertureInfo.planMU;


beamMU = [apertureInfo.beam(:).beamMU]';
beamArea = [apertureInfo.beam(:).beamArea]';
beamModulation = [apertureInfo.beam(:).beamModulation]';
beamK = [apertureInfo.beam(:).beamK]';


apertureMU(isnan(apertureMU)) = [];
apertureArea(isnan(apertureArea)) = [];


figure
histogram2(beamMU,beamArea,5,'DisplayStyle','tile')
xlabel('Monitor units')
ylabel('Beam area (cm^2)')
colorbar

figure
plot(apertureMU,apertureArea,'.')
xlabel('Monitor units')
ylabel('Aperture area (cm^2)')

figure
histogram2(beamMU,beamModulation,5,'DisplayStyle','tile')
xlabel('Monitor units')
ylabel('Beam modulation')
colorbar

figure
histogram2(beamArea,beamModulation,5,'DisplayStyle','tile')
xlabel('Area (cm^2)')
ylabel('Beam modulation')
colorbar

figure
plot(beamK,'.')
xlabel('Beam number')
ylabel('k = N*(1-BM)')

l = 0;
if pln.VMAT
    %All of these are vectors
    %Each entry corresponds to a beam angle
    %Later, we will convert these to histograms, find max, mean, min, etc.
    gantryRot = zeros(1,size(pln.optGantryAngles,2)-1);
    MURate = gantryRot;
    maxLeafSpeed = 0*pln.optGantryAngles;
    
    
    totTime = 0;
    
    for i = 1:size(apertureInfo.beam,2)
        if apertureInfo.beam(i).numOfShapes && l < apertureInfo.totalNumOfShapes-1 %only optimized beams have their time in the data struct
            %also, do not do the very last optimized beam!
            l = l+1;
            totTime = totTime+apertureInfo.beam(i).time; %time until next optimized beam
            gantryRot(l) = apertureInfo.beam(i).gantryRot;
            MURate(l) = apertureInfo.beam(i).MURate*60;
            maxLeafSpeed(l) = apertureInfo.beam(i).maxLeafSpeed/10;
        end
    end
    apertureInfo.time = totTime;
    
    %histogram of e.g., Modulation index vs gantry rotation, max leaf speed,
    %etc?
    figure
    hist(gantryRot)
    xlabel('Gantry rotation speed (deg/s)')
    figure
    hist(MURate)
    xlabel('Dose rate (MU/min)')
    figure
    hist(maxLeafSpeed)
    xlabel('Maximum leaf speed (cm/s)')
    
    fprintf('Min MU rate = %.1f',min(MURate));
end

result.apertureInfo = apertureInfo;

