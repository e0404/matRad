function matRad_calcVoxelWeighting(i,j,cst,d_i,d_ref,d_ref2)

global matRad_backprojectionFlag;
global matRad_voxelWeighting;

if matRad_backprojectionFlag 
    if matRad_voxelWeighting{i,2}
        
        if isequal(cst{i,6}(j).type,'min DCH objective') && d_ref < d_ref2 ||...
           isequal(cst{i,6}(j).type,'max DCH objective') && d_ref > d_ref2

        matRad_voxelWeighting{i,1} = 1;

        else
            
        % decimal place for dose comparison
        doseDecimalPlace = 1;        
       
        matRad_voxelWeighting{i,1} = zeros(1,numel(d_i));
            
        % apply lower and upper dose limits
        if isequal(cst{i,6}(j).type, 'max DCH objective')
            d_lower = d_ref;
            d_upper = d_ref2;
        elseif isequal(cst{i,6}(j).type, 'min DCH objective')
            d_lower = d_ref2;
            d_upper = d_ref;
        end
        
        logicalDoseLimits = (d_i >= round(d_lower*10^doseDecimalPlace)/(10^doseDecimalPlace) - 0.5*10^(-doseDecimalPlace) &...
                             d_i <= round(d_upper*10^doseDecimalPlace)/(10^doseDecimalPlace) + 0.4*10^(-doseDecimalPlace) );
                         
        if sum(logicalDoseLimits) > 0                 
            d_i               = d_i(logicalDoseLimits);
            minDistToVOI      = cst{i,5}.minDistToVOI(logicalDoseLimits);  

            % round dose values
            d_i = round(d_i*10^doseDecimalPlace)/(10^doseDecimalPlace);     

            % create logical dose mask (combine all dose entries)
            logicalDoseMask = bsxfun(@eq,sparse(d_i'),d_i);

            % apply mask on distances
            diagDist   = spdiags(minDistToVOI',0,numel(minDistToVOI),numel(minDistToVOI));
            eqDoseDist = diagDist*logicalDoseMask;

            % find min distances for every dose
            [~,jj] = find(eqDoseDist);
            weighting = accumarray(jj,nonzeros(eqDoseDist),[],@min);

            matRad_voxelWeighting{i,1}(logicalDoseLimits) = 1 + 4 * (weighting'./minDistToVOI);
        end
            
        end
        
    end
end


end
 