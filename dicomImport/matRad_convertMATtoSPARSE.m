function [ dij ] = matRad_convertMATtoSPARSE( dijFileNames,LETFileNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NumBeam         = length(dijFileNames);
NumEntries      = 0;

vDimesions.minX = inf;
vDimesions.maxX = 0;
vDimesions.minY = inf;
vDimesions.maxY = 0;
vDimesions.minZ = inf;
vDimesions.maxZ = 0;

cDim = {'X','Y','Z'};

for i = 1:NumBeam    
    
    load(dijFileNames{i});
    dijStruct.(['beam' num2str(i)]) = influenceStruct;    
    
    NumEntries = NumEntries + numel(influenceStruct.ixBeam);
    
    % check dimensions of sub cube
    for k = 1:numel(cDim)
        currMinField = (['min' cDim{1,k}]);
        currMaxField = (['max' cDim{1,k}]);
        if influenceStruct.stats.(currMinField) < vDimesions.(currMinField)
            vDimesions.(currMinField) = influenceStruct.stats.(currMinField);
        end
        if influenceStruct.stats.(currMaxField) > vDimesions.(currMaxField)
             vDimesions.(currMaxField) = influenceStruct.stats.(currMaxField);
        end
    end
    clear 'influenceStruct'   
end
 
newCubeDimension = [vDimesions.maxX-vDimesions.minX+1 vDimesions.maxY-vDimesions.minY+1 vDimesions.maxZ-vDimesions.minZ+1];
vOffset          = [vDimesions.minX vDimesions.minY vDimesions.minZ];

linIdx    = zeros(NumEntries,1);
beamletIx = zeros(NumEntries,1);
data      = zeros(NumEntries,1);
Cnt       = 1;
CntSpot   = 1;

for i = 1:NumBeam 
   NumEntriesCurrBeam = numel(dijStruct.(['beam' num2str(i)]).ixBeam);
   mCoord             =  bsxfun(@minus,dijStruct.(['beam' num2str(i)]).mSubScript,vOffset-1);
   linIdx(Cnt:Cnt+NumEntriesCurrBeam-1)    = sub2ind(newCubeDimension,mCoord(:,1),mCoord(:,2),mCoord(:,3));
   beamletIx(Cnt:Cnt+NumEntriesCurrBeam-1) = dijStruct.(['beam' num2str(i)]).ixBeamlet+CntSpot;
   data(Cnt:Cnt+NumEntriesCurrBeam-1)      = dijStruct.(['beam' num2str(i)]).data;
   CntSpot = CntSpot + min(dijStruct.(['beam' num2str(1)]).ixBeamlet);
   
   % missmatch of beamlet indices 
   aa = unique((dijStruct.(['beam' num2str(1)]).ixBeamlet));  
   bb = unique((dijStruct.(['beam' num2str(2)]).ixBeamlet));
   
   Cnt     = Cnt + NumEntriesCurrBeam;
end

clear mCoord dijStruct

dij.dose{1} = sparse(linIdx,beamletIx,data,prod(newCubeDimension),max(beamletIx));

%S = sparse(i,j,v,m,n)

% get dimension of sub cube ( dose cube




end

