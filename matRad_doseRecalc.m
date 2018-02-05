function recalc = matRad_doseRecalc(cst,pln,recalc,ct,apertureInfo)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to recalculate dose on a Dij angular resolution, or using
% the dynamic method, whichever the user wants%
%
% call
%   resultGUI = matRad_doseRecalc(dij,apertureInfo,resultGUI,pln)
%
% input
%   dij:            matRad dij struct
%   apertureInfo:   aperture shape info struct
%   resultGUI:      resultGUI struct to which the output data will be added, if
%                   this field is empty optResult struct will be created
%                   (optional)
%
% output
%   resultGUI:  struct containing optimized fluence vector, dose, and
%               shape info
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

recalc.apertureInfo = apertureInfo;
recalc.apertureInfo.recalcDynamic = recalc.dynamic;


%recalculate dose with finer gantry angles
%to do this, we need Dij matrices at these new angles
%Calculate dose directly

%first, we need to update/generate new apertures for these angles
recalc.pln = matRad_VMATGantryAngles(recalc.pln,'new');
if ~recalc.interpNew
    %we will not interpolate new apertures
    %easiest way to to this is to make ALL gantryAngles optGantryAngles
    recalc.pln.optGantryAngles = recalc.pln.gantryAngles;
end

cd stf
fname = sprintf('%.1f deg.mat',recalc.pln.minGantryAngleRes);
if exist(fname,'file')
    load(fname);
else
    stf = matRad_generateStf(ct,cst,recalc.pln);
    save(fname,'stf')
end
recalc.stf = stf;
clear stf
cd ..

if ~recalc.interpNew || ~recalc.dijNew
    %duplicate any beam angles that are directly between two old
    %ones
    duplicate = false(size(recalc.pln.gantryAngles));
    for i = 1:numel(recalc.pln.gantryAngles)
        if numel(find(abs(recalc.pln.gantryAngles(i)-pln.gantryAngles) == min(abs(recalc.pln.gantryAngles(i)-pln.gantryAngles)))) > 1
            duplicate(i) = true;
        end
    end
    newGantryAngles = zeros(1,numel(recalc.pln.gantryAngles)+nnz(duplicate));
    newCouchAngles = zeros(1,numel(recalc.pln.gantryAngles)+nnz(duplicate));
    tempStf = recalc.stf;
    recalc.stf(1).copyInd = [];
    tempStf(1).copyInd = [];
    recalc.stf(1).stfCorr = [];
    tempStf(1).stfCorr = [];
    j = 1;
    for i = 1:numel(recalc.pln.gantryAngles)
        if duplicate(i)
            tempStf(j).stfCorr = false;
            newGantryAngles(j) = recalc.pln.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.couchAngles(i);
            tempStf(j) = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i-1).gantryAngle;
            tempStf(j).copyInd = 1;
            
            j = j+1;
            
            newGantryAngles(j) = recalc.pln.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.couchAngles(i);
            tempStf(j) = recalc.stf(i);
            tempStf(j).gantryAngle = recalc.stf(i+1).gantryAngle;
            tempStf(j).copyInd = 2;
        else
            tempStf(j).stfCorr = true;
            newGantryAngles(j) = recalc.pln.gantryAngles(i);
            newCouchAngles(j) = recalc.pln.couchAngles(i);
            tempStf(j) = recalc.stf(i);
        end
        j = j+1;
    end
    recalc.pln.gantryAngles = newGantryAngles;
    recalc.pln.couchAngles = newCouchAngles;
    recalc.pln.numOfBeams = numel(recalc.pln.gantryAngles);
    %recalc.pln.optGantryAngles = recalc.pln.gantryAngles;
    recalc.stf = tempStf;
end

recalc = matRad_recalcApertureInfo(recalc,recalc.apertureInfo);

if ~recalc.interpNew || ~recalc.dijNew
    tempPln = recalc.pln;
    tempStf = recalc.stf;
    for i = 1:numel(tempPln.gantryAngles)
        diff = abs(tempPln.gantryAngles(i)-pln.gantryAngles);
        minDiffInd = diff == min(diff);
        minDiffInd1 = find(tempPln.gantryAngles == pln.gantryAngles(find(minDiffInd,1,'first')));
        minDiffInd2 = find(tempPln.gantryAngles == pln.gantryAngles(find(minDiffInd,1,'last')));
        
        if ~recalc.dijNew
            if isempty(recalc.stf(i).copyInd)
                recalc.stf(i) = tempStf(minDiffInd1);
                recalc.pln.gantryAngles(i) = tempPln.gantryAngles(minDiffInd1);
            elseif recalc.stf(i).copyInd == 1
                recalc.stf(i) = tempStf(minDiffInd1);
                recalc.pln.gantryAngles(i) = tempPln.gantryAngles(minDiffInd1);
            elseif recalc.stf(i).copyInd == 2
                recalc.stf(i) = tempStf(minDiffInd2);
                recalc.pln.gantryAngles(i) = tempPln.gantryAngles(minDiffInd2);
            end
        elseif ~recalc.interpNew
            if numel(minDiffInd) > 1
                recalc.stf(i).gantryAngle = tempPln.gantryAngles(i);
            end
        end
    end
end

if recalc.dynamic
    recalc.apertureInfo =  matRad_daoVec2ApertureInfo_VMATrecalcDynamic(recalc.apertureInfo,recalc.apertureInfo.apertureVector);
else
    recalc.apertureInfo =  matRad_daoVec2ApertureInfo_VMATstatic(recalc.apertureInfo,recalc.apertureInfo.apertureVector);
end

recalc.resultGUI = matRad_calcDoseDirect(ct,recalc.stf,recalc.pln,cst,recalc.apertureInfo.bixelWeights);








