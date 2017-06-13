sigma_ray = sigma_rayvec(5);
for i=1:30
    radiusmat(:,i) = linspace(sigma_ray/4, 0.95*sigma_ray/1.1, 20);
end
for i=1:20
    sigma_submat(i,:) = linspace(1.1*radiusmat(i,1),1.5*sigma_ray,30);
end

H = reshape(maxi_rad(:,:,5),[1,600]);
[Hsort,Idx] = sort(H,2);
sigma_submat_sort = sigma_submat(Idx);
radiusmat_sort = radiusmat(Idx);
X_rad_temp = X_rad(:,:,5,2); sigma_res = X_rad_temp(Idx);
X_rad_temp = X_rad(:,:,5,1); Ampl_res = X_rad_temp(Idx);

res_idx = sigma_res>=sqrt(7^2-sigma_submat_sort.^2) & sigma_submat_sort<=0.95*sigma_ray ;
superidx = find(res_idx,1);

Untitled2(sigma_ray,sigma_submat_sort(superidx),radiusmat_sort(superidx),n,...
    [Ampl_res(superidx) sigma_res(superidx)], 'circle')
%%
m=0;
for i=1:12
    for j=1:10
        if sigma_subvec(j) >= radiusvec(i)
            m=m+1;
            coord(1,m) =  sigma_subvec(j);
            coord(2,m) =  radiusvec(i);
        end
    end
end

subplot(1,2,1)
surf(sigma_subvec,radiusvec,maxi_rad)
title('max error')
xlabel('\sigma_{sub}'); ylabel('radius')
view([0 90])

subplot(1,2,2)
surf(sigma_subvec,radiusvec,X_rad(:,:,:,2))
title('Sigma')
xlabel('\sigma_{sub}'); ylabel('radius')
view([0 90])

hold
scatter3(coord(1,:),coord(2,:),10.*ones([1 length(coord(1,:))])','xr','LineWidth',10)











