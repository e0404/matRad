function cst = matRad_addVOIRing(cst,ct,myMargin,calcDistFlag)
 
Counter = 0;

for  i = 1:size(cst,1)
    if sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0
        
        Counter = Counter + 1;  

        % generate VOI cube
        V          = cst{i,4}{1};
        VOICube    = zeros(ct.cubeDim);
        VOICube(V) = 1;

        % add VOI ring
        VOICubeWithRing = matRad_addMargin(VOICube,cst,ct.resolution,myMargin,true);
        VwithRing       = find(VOICubeWithRing>0);

        % create cst with ring structure
        cstRing{Counter,1}    = size(cst,1) - 1 + Counter;
        cstRing{Counter,2}    = [cst{i,2},'Ring'];
        cstRing{Counter,3}    = cst{i,3};
        cstRing{Counter,4}{1} = setdiff(VwithRing,V);
        cstRing{Counter,5}    = cst{i,5};

        % calc min distance to VOI for every VOI ring voxel
        if calcDistFlag
            cstRing{Counter,5}.minDistToTarget = matRad_calcMinDist(ct,cstRing{Counter,4}{1},cst{i,4}{1});
        end

        % pass coverage based objective/constraint specification to VOI ring structure
        logidx             = strcmp({cst{i,6}(:).robustness},'coverage');
        cstRing{Counter,6} = cst{i,6}(logidx);
        cst{i,6}           = cst{i,6}(~logidx);
    end  
end

cst = [cst;cstRing];

end