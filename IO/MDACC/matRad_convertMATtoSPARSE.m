function [ dij ] = matRad_convertMATtoSPARSE(ijFileNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


NumBeam         = cell2mat(cellfun(@(x) sum(ismember(ijFileNames(:,2),x)),unique(ijFileNames(:,2),'stable'),'un',0));
NumEntries      = 0;

vDimesions.minX = inf;
vDimesions.maxX = 0;
vDimesions.minY = inf;
vDimesions.maxY = 0;
vDimesions.minZ = inf;
vDimesions.maxZ = 0;

cDim = {'X','Y','Z'};


for j = 1:size(ijFileNames,1)

    for i = 1:NumBeam    

        load(ijFileNames{i,1});
        dijStruct.(['beam' num2str(i)]) = influenceStruct;    

        NumEntries = NumEntries + numel(influenceStruct.ixBeam);

        % check extension of each dimension of the sub dose cube
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

    linIdx    = zeros(NumEntries,1,'double');
    beamletIdx = zeros(NumEntries,1,'double');
    data      = zeros(NumEntries,1,'double');
    
    Cnt       = 1;
    CntSpot   = 1;

    for i = 1:NumBeam 
       NumEntriesCurrBeam                      = numel(dijStruct.(['beam' num2str(i)]).ixBeam);
       
       mCoord                                  =  bsxfun(@minus,dijStruct.(['beam' num2str(i)]).mSubScript,vOffset-1);
       dijStruct.(['beam' num2str(i)])         = rmfield(dijStruct.(['beam' num2str(i)]),'mSubScript');
       
       linIdx(Cnt:Cnt+NumEntriesCurrBeam-1)    = sub2ind(newCubeDimension,mCoord(:,1),mCoord(:,2),mCoord(:,3));
       
       beamletIdx(Cnt:Cnt+NumEntriesCurrBeam-1) = dijStruct.(['beam' num2str(i)]).ixBeamlet+double(CntSpot);
       dijStruct.(['beam' num2str(i)])         = rmfield(dijStruct.(['beam' num2str(i)]),'ixBeamlet');
       
       data(Cnt:Cnt+NumEntriesCurrBeam-1)      = dijStruct.(['beam' num2str(i)]).data;
       dijStruct.(['beam' num2str(i)])         = rmfield(dijStruct.(['beam' num2str(i)]),'data');
       
       CntSpot = CntSpot + max(beamletIdx(Cnt:Cnt+NumEntriesCurrBeam-1));

       Cnt     = Cnt + NumEntriesCurrBeam;
    end

    clear mCoord dijStruct
    
    if strcmp(ijFileNames(j,2),'Dij')
        dij.dose{1} = sparse(linIdx,beamletIdx,data,prod(newCubeDimension),max(beamletIdx));
    elseif strcmp(ijFileNames{:,2},'LETij')
        dij.LET{1} = sparse(linIdx,beamletIdx,data,prod(newCubeDimension),max(beamletIdx));
    end

    clear 'linIdx' 'beamletIdx' 'data'

end
%S = sparse(i,j,v,m,n)

% get dimension of sub cube ( dose cube




end

