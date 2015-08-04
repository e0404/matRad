function [Ix_l, Ix_r] = tk_getLeafPos(shapeMap)

[dimZ, dimX] = size(shapeMap);
Ix_l = NaN * ones(dimZ,1);
Ix_r = NaN * ones(dimZ,1);

% loop over all lines of the Map
for k=1:size(shapeMap,1)
    % get one line at a time
    tempLine = shapeMap(k,:);
    
    %% search for the coordinates of the right (Ix_r) and left (Ix_l) leaf
        bixelSum = sum(tempLine(:));
    
        % 1. check if any bixel is open
        if  bixelSum == 0
            Ix_l(k) = round(dimX/2);
            Ix_r(k) = round(dimX/2);
            continue
        end
        
        % 2. check if all bixels are open
        if bixelSum == dimX
            Ix_l(k) = 0;
            Ix_r(k) = dimX;
            continue
        end            
        
        % 3. find out individual leaf positions
            % check if left leaf is completely open
            if tempLine(1) == 1
                Ix_l(k) = 0;
                Ix_r(k) = bixelSum;
                continue
            end
            % check if right leaf is completely open
            if tempLine(dimX) == 1
                Ix_l(k) = dimX - bixelSum;
                Ix_r(k) = dimX;
                continue
            end
            
            % check for firt open bixel
            for l=2:dimX-1
                if tempLine(l) == 1
                    Ix_l(k) = l-1;
                    Ix_r(k) = l-1 + dimX;
                    continue
                end           
            end               
    
end

end