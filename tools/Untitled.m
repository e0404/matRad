%surf(radiusvec,sigma_subvec,reshape(maxi_rad(3,:,:),[20 20]))

for i=4
    sigma_ray = sigma_rayvec(i);
    radiusvec = linspace(sigma_ray/4, 0.95*sigma_ray/1.1, 20);
    minerr(i) = min(min(maxi_rad(:,:,i)));
    [alpha,beta] = find(maxi_rad(:,:,i)==minerr(i));
    minrad(i) = radiusvec(alpha);
    sigma_subvec = linspace(1.1*minrad(i),1.5*sigma_ray,30);
    sigma_subvec = sort(sigma_subvec);
    minsig(i) = sigma_subvec(beta);
    disp([minrad(i) minsig(i)])
    figure(1)
    hold off
    surf(sigma_subvec,radiusvec,maxi_rad(:,:,i))
    hold
    scatter(sigma_subvec(beta),radiusvec(alpha),'xr','LineWidth',10)
    x = [0.95*sigma_ray 0.95*sigma_ray];
    y = [sigma_ray/4 0.95*sigma_ray/1.1];
    z = [0 40; 0 40];
    surf(x,y,z,'LineWidth',.1,'EdgeAlpha',.3,'FaceAlpha',.3,'FaceColor','g')
    sigma_subvec2 = sigma_subvec(sigma_subvec < 0.95*sigma_ray);
    minerrcorr(i) = min(min(maxi_rad(:,length(sigma_subvec2),i)));
    [alpha,beta] = find(maxi_rad(:,:,i)==minerrcorr(i));
    minradcorr(i) = radiusvec(alpha);
    sigma_subvec = linspace(1.1*minradcorr(i),1.5*sigma_ray,30);
    sigma_subvec2 = sigma_subvec(sigma_subvec < 0.95*sigma_ray);
    sigma_subvec2 = sort(sigma_subvec2);
    minsigcorr(i) = sigma_subvec2(beta);
    disp([i minradcorr(i) minsigcorr(i)])
    scatter(sigma_subvec(beta),radiusvec(alpha),'xr','LineWidth',10)
    figure
    %hold off
    scatter3(0,0,0)
    Untitled2(minradcorr(i),sigma_subvec2(beta),radiusvec(alpha),2,X_rad(alpha,beta,1,:),'circle')
    Xcorr(i,:) = X_rad(alpha,beta,1,:);
%     figure(2)
%     hold off
%     X_rad(alpha,beta,1,:)
    
end