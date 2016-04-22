function cst = matRad_addVOIRing(cst,ct,myMargin,calcWeightingFlag,voxelWeightingType,multScen)

global matRad_voxelWeighting;
 
covFlag = 0;
Counter = 0;

for  i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        if sum(strcmp({cst{i,6}(:).robustness},'coverage')) > 0

            covFlag = 1;
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
            if calcWeightingFlag
                if isequal(voxelWeightingType,'heurWeighting')
                    
                    [cstRing{Counter,5}.minDistToVOI,~] = matRad_calcMinDist(ct,cstRing{Counter,4}{1},cst{i,4}{1});
                    matRad_voxelWeighting{Counter + size(cst,1),1}  = 1;
                    matRad_voxelWeighting{Counter + size(cst,1),2}  = true;
                    
                elseif isequal(voxelWeightingType,'probWeighting')
                    
                    [~,minDistToVOI] = matRad_calcMinDist(ct,cstRing{Counter,4}{1},cst{i,4}{1});
                    sigma            = multScen.shiftSize';
                    
                    if sum(sigma > 0) == 0
                        error('sigma_x, sigma_y and sigma_z equals to zero, prob voxel weighting impossible')
                        
                    elseif sum(sigma > 0) == 1
                        normalization = 1/(sigma(sigma > 0)*sqrt(2*pi));
                        
                    elseif sum(sigma > 0) == 2
                        normalization = 1/(prod(sigma(sigma > 0))*2*pi); 
                        
                    elseif sum(sigma > 0) == 3
                        normalization = 1/(prod(sigma)*(2*pi)^(3/2));
                        
                    end

                    cstRing{Counter,5}.voxelWeighting = 5*normalization* exp(-sum(minDistToVOI(sigma > 0,:).^2./repmat(2*sigma(sigma > 0).^2,1,size(minDistToVOI,2)),1));
                end
            end

            % pass coverage based objective/constraint specification to VOI ring structure
            logidx             = strcmp({cst{i,6}(:).robustness},'coverage');
            cstRing{Counter,6} = cst{i,6}(logidx);
            cst{i,6}           = cst{i,6}(~logidx);
        end 
    end
end

if covFlag
    cst = [cst;cstRing];
end

end