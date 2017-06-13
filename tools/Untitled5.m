

n=2;
method='circle';
sigma_rayvec = 6;
figure
i=0;

for s=1:1
    sigma_ray = sigma_rayvec(s);
    radiusvec = linspace(sigma_ray/2, 0.95*sigma_ray/1.1, 12);
    for k=1:12
        radius = radiusvec(k);
        sigma_subvec = linspace(1.1*sigma_ray/2,0.95*sigma_ray,10);
        for m=1:10
            sigma_sub = sigma_subvec(m);
            
            X1 = X_rad(k,m,s,:);
            
            if X1(2) >= radius && maxi_rad(k,m) <= 1
                i=i+1;
                %Untitled2(sigma_ray,sigma_sub,radius,n,X1,method,maxi_rad(k,m))
                tripletta(i,:) = ([maxi_rad(k,m), radius, sigma_sub, X1(1), X1(2)]);
            end
        end
    end
end

tripletta(tripletta(:,1)==min(tripletta(:,1)),:)
