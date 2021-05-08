function posAtSlice = getNeedlePosition(obj,rpItem,posZ,noNeedle)
%GETNEEDLEPOSITION Summary of this function goes here
%   Detailed explanation goes here

prec = 1e-9;


if nargin>3 && isempty(posZ)
    posAtSlice = obj.RP.(rpItem).seeds.dwellPosition{noNeedle};
elseif nargin<4
    posAtSlice = zeros(size(obj.RP.(rpItem).seeds.dwellPosition,1),3);
    for in = 1:size(obj.RP.(rpItem).seeds.dwellPosition,1)
        dwellPos = obj.RP.(rpItem).seeds.dwellPosition{in};
        
        if posZ<min(dwellPos(:,3))
            dz = abs(dwellPos(:,3)-min(dwellPos(:,3)));
            ind = dz<prec;
            posAtSlice(1,:) = [dwellPos(ind,1:2),posZ];
        elseif posZ>max(dwellPos(:,3))
            dz = abs(dwellPos(:,3)-max(dwellPos(:,3)));
            ind = dz<prec;
            posAtSlice(1,:) = [dwellPos(ind,1:2),posZ];
        else
            
            dz = abs(dwellPos(:,3)-posZ);
            ind = dz<prec;
            if any(ind)
                posAtSlice(in,:)   = dwellPos(ind,:);
            else
                [~,indZ] = sort(dz);
                
                p1 = dwellPos(indZ(1),:);
                p2 = dwellPos(indZ(2),:);
                
                lambda             = (posZ-p1(3))/(p2(3)-p1(3));
                posAtSlice(in,:)   = p1 + lambda*(p2-p1);
            end
        end
    end
else
    if ~isempty(obj.RP.(rpItem).seeds.dwellPosition{noNeedle})
        posAtSlice = zeros(length(posZ),3);
        dwellPos = obj.RP.(rpItem).seeds.dwellPosition{noNeedle};
        for iz=1:length(posZ)
            
            if posZ(iz)<min(dwellPos(:,3))
                dz = abs(dwellPos(:,3)-min(dwellPos(:,3)));
                ind = dz<prec;
                posAtSlice(iz,:) = [dwellPos(ind,1:2),posZ(iz)];
            elseif posZ(iz)>max(dwellPos(:,3))
                dz = abs(dwellPos(:,3)-max(dwellPos(:,3)));
                ind = dz<prec;
                posAtSlice(iz,:) = [dwellPos(ind,1:2),posZ(iz)];
            else
                dz = abs(dwellPos(:,3)-posZ(iz));
                ind = dz<prec;
                if any(ind)
                    posAtSlice(iz,:) = dwellPos(ind,:);
                else
                    [~,indZ] = sort(dz);
                    
                    p1 = dwellPos(indZ(1),:);
                    p2 = dwellPos(indZ(2),:);
                    
                    lambda     = (posZ(iz)-p1(3))/(p2(3)-p1(3));
                    posAtSlice(iz,:) = p1 + lambda*(p2-p1);
                end
            end
        end
    else
        posAtSlice = [];
    end
end
end


