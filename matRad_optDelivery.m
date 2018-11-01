function result = matRad_optDelivery(result,fast)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad: optimize VMAT delivery
% 
% call
%   matRad_optDelivery(result,fast)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   fast:               1 => fastest possible delivery
%                       0 => mutliply delivery time by 10%
%
% output
%   apertureInfo:       aperture shape info struct
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

if nargin < 3
    fast = 1;
end


%speed up delivery time, when it is permitted by constraints
%constraints to consider: doserate, leaf speed, and gantry speed

%Do this after DAO

apertureInfo = result.apertureInfo;

fileName = apertureInfo.propVMAT.machineConstraintFile;
try
    load(fileName,'machine');
catch
    error(['Could not find the following machine file: ' fileName ]);
end



%calculate max leaf speed
apertureInfo = matRad_maxLeafSpeed(apertureInfo);

doInterp = 0;

for i = 1:size(apertureInfo.beam,2)
    if apertureInfo.propVMAT.beam(i).DAOBeam
        
        %all of these should be greater than 1, since DAO respects the
        %constraints
        
        %if one of them is less than 1, then a constraint is violated
        factorMURate    = machine.constraints.monitorUnitRate(2)/apertureInfo.beam(i).shape(1).MURate;
        factorLeafSpeed = machine.constraints.leafSpeed(2)/apertureInfo.beam(i).maxLeafSpeed;
        factorGantryRot = machine.constraints.gantryRotationSpeed(2)/apertureInfo.beam(i).gantryRot;
        
        %The constraint that is limiting the speed the most is the one
        %whose factor is closest to 1
        factor = min([factorMURate factorLeafSpeed factorGantryRot]);
        if ~fast
            %if the limiting rate is already 10% lower than the limit,
            %then do nothing (factor = 1)
            %otherwise, scale rates so that the limiting rate is 10% lower
            %than the limit
            factor = min([1 factor*0.9]);
        end
        
        %multiply each speed by this factor
        apertureInfo.beam(i).shape.MURate = factor*apertureInfo.beam(i).shape.MURate;
        apertureInfo.beam(i).maxLeafSpeed = factor*apertureInfo.beam(i).maxLeafSpeed;
        apertureInfo.beam(i).gantryRot = factor*apertureInfo.beam(i).gantryRot;
        apertureInfo.beam(i).time = apertureInfo.beam(i).time/factor;
        
        factorMURate = machine.constraints.monitorUnitRate(1)/apertureInfo.beam(i).shape(1).MURate;
        
        if factorMURate > 1
            apertureInfo.beam(i).shape(1).MURate = factorMURate*apertureInfo.beam(i).shape(1).MURate;
            apertureInfo.beam(i).shape(1).weight = factorMURate*apertureInfo.beam(i).shape(1).weight;
            
            doInterp = 1;
        end
    end
end

%recalculate vector with new times

%%%LOOK PAST THIS
[apertureInfo.apertureVector,~,~] = matRad_daoApertureInfo2Vec(apertureInfo);

%redo interpolation
apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfo.apertureVector);

if doInterp
    fprintf('\n\nWE ARE REDOING INTERPOLATION\n\n');
else
    fprintf('\n\nWE ARE NOT REDOING INTERPOLATION\n\n');
end

result.apertureInfo = apertureInfo;

