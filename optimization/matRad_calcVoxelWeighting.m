function matRad_calcVoxelWeighting(i,j,cst,d_i,d_ref,d_ref2)

global matRad_backprojectionFlag;
global matRad_iteration;
global matRad_voxelWeighting;

if matRad_backprojectionFlag 
    if matRad_voxelWeighting{i,2}
        
        if isequal(cst{i,6}(j).type,'min DCH objective') && d_ref < d_ref2 ||...
           isequal(cst{i,6}(j).type,'max DCH objective') && d_ref > d_ref2 ||...
           matRad_iteration < 5

        matRad_voxelWeighting{i} = 1;

        else
        
        % round dose values
        d_i = round(d_i*10)/10;     
            
        % create logical dose mask (combine all dose entries)
        logicalDoseMask = bsxfun(@eq,sparse(d_i'),d_i);

        % apply mask on distances
        diagDist   = spdiags(cst{i,5}.minDistToVOI',0,numel(cst{i,5}.minDistToVOI),numel(cst{i,5}.minDistToVOI));
        eqDoseDist = diagDist*logicalDoseMask;

        % find min distances for every dose
        [~,jj] = find(eqDoseDist);
        weighting = accumarray(jj,nonzeros(eqDoseDist),[],@min);

        matRad_voxelWeighting{i} = 1 + 4 * (weighting'./cst{i,5}.minDistToVOI);

%         for k = 1:length(d_i)
%             weighting(k) = min(cst{i,5}.minDistToVOI(d_i(k) == d_i));
%         end
%         matRad_voxelWeighting{i} = 1 + 4 * (weighting./cst{i,5}.minDistToVOI);
            
        end
        
    end
end


end
 