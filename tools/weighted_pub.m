%% mR_findWeights.m
% This function provides the weights of sub-beams for different standard
% deviations of the lateral spread of the incoming beam

%%
% In the plots we can see: the percentage difference between the gaussian
% representig the lateral spread and the weighted sum of gaussians.
% sigma_ray = 7;
% n = 2;
% radius = 4;
% sigma_sub = 5.5;
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'square');
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'circle');

%%
% sigma_ray = 5;
% n = 2;
% radius = 3;
% sigma_sub = 4;
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'square');
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'circle');
% 
% %%
% sigma_ray = 9;
% n = 3;
% radius = 5;
% sigma_sub = 7;
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'square');
% X1 = mR_findWeights(sigma_ray, sigma_sub, n, radius, 'circle');

%% Circular weights
% We have here an analysis of the weighted sum using a circular shape.
% 
load('pezzWeightsData_circle') % gotta save to it and change name
method = 'circle';
clear tripletta
clear result

for s=1:13
    i=0;
    sigma_ray = sigma_rayvec(s);
    radiusvec = linspace(sigma_ray/2, 0.95*sigma_ray/1.1, 22);
    for k=1:22
        radius = radiusvec(k);
        sigma_subvec = linspace(1.1*sigma_ray/2,0.95*sigma_ray,20);
        for m=1:20
            sigma_sub = sigma_subvec(m);
            
            X1 = X_rad(k,m,s,:);
            
            if X1(2) >= radius && maxi_rad(k,m,s) <= 1 % this comes from 
                % the sampling theorem
                i=i+1;
                %Untitled2(sigma_ray,sigma_sub,radius,n,X1,method,maxi_rad(k,m))
                tripletta{s}(i,:) = ([maxi_rad(k,m,s), radius, sigma_sub, X1(1), X1(2), timei_rad(k,m,s)]);
            end
        end
    end
    result(s,:) = tripletta{s}(tripletta{s}(:,1)==min(tripletta{s}(:,1)),:);
%      disp(s)
%     disp(result)
%     Untitled2(sigma_ray, result(s,3), result(s,2), n, [result(s,4) result(s,5)], method, result(s,1))
end
s0 = find(result(:,1)==min(result(:,1)));
Untitled2(sigma_rayvec(s0), result(s0,3), result(s0,2), n, [result(s0,4) result(s0,5)], method, result(s0,1))

%%
% Here we plot the minimum values for each sigma_ray and we analyize some
% more info as time spent and other parameters
x = 2:0.01:16;

figure
plot(sigma_rayvec, result(:,1))
title('\sigma_{ray} Vs. Maximum (%) error')

figure
scatter(sigma_rayvec, result(:,2))
hold
coeff_rad_circ = polyfit(sigma_rayvec(1:10)',result(1:10,2),1);
plot(x, coeff_rad_circ(1).*x+coeff_rad_circ(2),'r')
title('\sigma_{ray} Vs. radius (mm)')

figure
scatter(sigma_rayvec, result(:,3))
hold
coeff_sig_circ = polyfit(sigma_rayvec(1:10)',result(1:10,3),1);
plot(x, coeff_sig_circ(1).*x+coeff_sig_circ(2),'r')
title('\sigma_{ray} Vs. \sigma_{sub}')

figure
plot(sigma_rayvec, result(:,6))
title('\sigma_{ray} Vs. computation time')

figure
scatter(sigma_rayvec, result(:,4))
hold
coeff= polyfit(sigma_rayvec(1:10)',result(1:10,4),2);
plot(x, coeff(1).*x.^2+coeff(2).*x+coeff(3),'r')
title('\sigma_{ray} Vs. normalization of weights')

figure
scatter(sigma_rayvec, result(:,5))
hold
coeff= polyfit(sigma_rayvec(1:10)',result(1:10,5),1);
plot(x, coeff(1).*x+coeff(2),'r')
title('\sigma_{ray} Vs. \sigma_w')

coeff_rad_circ_correct = polyfit(sigma_rayvec(1:10)',result(1:10,2),1);
coeff_sig_circ_correct = polyfit(sigma_rayvec(1:10)',result(1:10,3),1);
coeff_sigW_circ_correct= polyfit(sigma_rayvec(1:10)',result(1:10,5),1);
coeff_w_circ_correct= polyfit(sigma_rayvec(1:10)',result(1:10,4),2);

%% Square shape
load('pezzWeightsData_square') % gotta save to it and change name
method = 'square';
clear tripletta
clear result

for s=1:13
    i=0;
    sigma_ray = sigma_rayvec(s);
    radiusvec = linspace(sigma_ray/2, 0.95*sigma_ray/1.1, 15);
    for k=1:15
        radius = radiusvec(k);
        sigma_subvec = linspace(1.1*sigma_ray/2,0.95*sigma_ray,12);
        for m=1:12
            sigma_sub = sigma_subvec(m);
            
            X1 = X_rad(k,m,s,:);
            
            if X1(2) >= radius && maxi_rad(k,m,s) <= 1
                i=i+1;
                %Untitled2(sigma_ray,sigma_sub,radius,n,X1,method,maxi_rad(k,m))
                tripletta{s}(i,:) = ([maxi_rad(k,m,s), radius, sigma_sub, X1(1), X1(2), timei_rad(k,m,s)]);
            end
        end
    end
    result(s,:) = tripletta{s}(tripletta{s}(:,1)==min(tripletta{s}(:,1)),:);
%     disp(s)
%     disp(result)
    %Untitled2(sigma_ray, result(s,3), result(s,2), n, [result(s,4) result(s,5)], method, result(s,1))
%     fprintf('\n')
%     %
end
s0 = find(result(:,1)==min(result(:,1)));
Untitled2(sigma_rayvec(s0), result(s0,3), result(s0,2), n, [result(s0,4) result(s0,5)], method, result(s0,1))

%%
figure
plot(sigma_rayvec, result(:,1))
title('\sigma_{ray} Vs. Maximum (%) error')

figure
scatter(sigma_rayvec, result(:,2))
hold
coeff_rad_sqr = polyfit(sigma_rayvec(1:12)',result(1:12,2),1);
plot(x, coeff_rad_sqr(1).*x+coeff_rad_sqr(2),'r')
title('\sigma_{ray} Vs. radius (mm)')

figure
scatter(sigma_rayvec, result(:,3))
hold
coeff_sig_sqr = polyfit(sigma_rayvec(1:12)',result(1:12,3),1);
plot(x, coeff_sig_sqr(1).*x+coeff_sig_sqr(2),'r')
title('\sigma_{ray} Vs. \sigma_{sub}')

figure
plot(sigma_rayvec, result(:,6))
title('\sigma_{ray} Vs. computation time')

figure
scatter(sigma_rayvec, result(:,4))
hold
coeff= polyfit(sigma_rayvec(1:12)',result(1:12,4),2);
plot(x, coeff(1).*x.^2+coeff(2).*x+coeff(3),'r')
title('\sigma_{ray} Vs. normalization of weights')

figure
scatter(sigma_rayvec, result(:,5))
hold
coeff= polyfit(sigma_rayvec(1:12)',result(1:12,5),1);
plot(x, coeff(1).*x+coeff(2),'r')
title('\sigma_{ray} Vs. \sigma_w')
% 
coeff_rad_sqr_correct = polyfit(sigma_rayvec(1:12)',result(1:12,2),1);
coeff_sig_sqr_correct = polyfit(sigma_rayvec(1:12)',result(1:12,3),1);
coeff_sigW_sqr_correct= polyfit(sigma_rayvec(1:12)',result(1:12,5),1);
coeff_w_sqr_correct= polyfit(sigma_rayvec(1:12)',result(1:12,4),2);

%% n=3

load('pezzWeightsData_circle3') 
method = 'circle';
clear tripletta
clear result
n=3;

for s=1:12
    i=0;
    sigma_ray = sigma_rayvec(s);
    radiusvec = linspace(sigma_ray/2, 0.95*sigma_ray/1.1, 22);
    for k=1:22
        radius = radiusvec(k);
        sigma_subvec = linspace(1.1*sigma_ray/2,0.95*sigma_ray,20);
        for m=1:20
            sigma_sub = sigma_subvec(m);
            
            X1 = X_rad(k,m,s,:);
            
            if X1(2) >= radius && maxi_rad(k,m,s) <= 5 % this comes from 
                % the sampling theorem
                i=i+1;
                %Untitled2(sigma_ray,sigma_sub,radius,n,X1,method,maxi_rad(k,m))
                tripletta{s}(i,:) = ([maxi_rad(k,m,s), radius, sigma_sub, X1(1), X1(2), timei_rad(k,m,s)]);
            end
        end
    end
    result(s,:) = tripletta{s}(tripletta{s}(:,1)==min(tripletta{s}(:,1)),:);
%       disp(s)
%     disp(result)
%     Untitled2(sigma_ray, result(s,3), result(s,2), n, [result(s,4) result(s,5)], method, result(s,1))
end
s0 = find(result(:,1)==min(result(:,1)));
Untitled2(sigma_rayvec(s0), result(s0,3), result(s0,2), n, [result(s0,4) result(s0,5)], method, result(s0,1))

%%
% Here we plot the minimum values for each sigma_ray and we analyize some
% more info as time spent and other parameters
x = 7:0.01:20;

figure
plot(sigma_rayvec(1:12), result(:,1))
title('\sigma_{ray} Vs. Maximum (%) error')

figure
scatter(sigma_rayvec(1:12), result(:,2))
hold
coeff_rad_circ = polyfit(sigma_rayvec(1:12)',result(:,2),1);
plot(x, coeff_rad_circ(1).*x+coeff_rad_circ(2),'r')
title('\sigma_{ray} Vs. radius (mm)')

figure
scatter(sigma_rayvec(1:12), result(:,3))
hold
coeff_sig_circ = polyfit(sigma_rayvec(1:10)',result(1:10,3),1);
plot(x, coeff_sig_circ(1).*x+coeff_sig_circ(2),'r')
title('\sigma_{ray} Vs. \sigma_{sub}')

figure
plot(sigma_rayvec(1:12), result(:,6))
title('\sigma_{ray} Vs. computation time')

figure
scatter(sigma_rayvec(1:12), result(:,4))
hold
coeff= polyfit(sigma_rayvec(1:12)',result(1:12,4),2);
plot(x, coeff(1).*x.^2+coeff(2).*x+coeff(3),'r')
title('\sigma_{ray} Vs. normalization of weights')

figure
scatter(sigma_rayvec(1:12), result(:,5))
hold
coeff= polyfit(sigma_rayvec(1:10)',result(1:10,5),1);
plot(x, coeff(1).*x+coeff(2),'r')
title('\sigma_{ray} Vs. \sigma_w')

coeff_rad_circ3_correct = polyfit(sigma_rayvec(1:12)',result(1:12,2),1);
coeff_sig_circ3_correct = polyfit(sigma_rayvec(1:10)',result(1:10,3),1);
coeff_sigW_circ3_correct= polyfit(sigma_rayvec(1:10)',result(1:10,5),1);
coeff_w_circ3_correct= polyfit(sigma_rayvec(1:12)',result(1:12,4),2);
%% Square shape n=3
load('pezzWeightsData_square3') % gotta save to it and change name
method = 'square';
clear tripletta
clear result
n=3;

for s=1:10
    i=0;
    sigma_ray = sigma_rayvec(s);
    radiusvec = linspace(sigma_ray/2, 0.95*sigma_ray/1.1, 15);
    for k=1:15
        radius = radiusvec(k);
        sigma_subvec = linspace(1.1*sigma_ray/2,0.95*sigma_ray,12);
        for m=1:12
            sigma_sub = sigma_subvec(m);
            
            X1 = X_rad(k,m,s,:);
            
            if X1(2) >= radius && maxi_rad(k,m,s) <= 1
                i=i+1;
                %Untitled2(sigma_ray,sigma_sub,radius,n,X1,method,maxi_rad(k,m))
                tripletta{s}(i,:) = ([maxi_rad(k,m,s), radius, sigma_sub, X1(1), X1(2), timei_rad(k,m,s)]);
            end
        end
    end
    result(s,:) = tripletta{s}(tripletta{s}(:,1)==min(tripletta{s}(:,1)),:);
%     disp(s)
%     disp(result)
    %Untitled2(sigma_ray, result(s,3), result(s,2), n, [result(s,4) result(s,5)], method, result(s,1))
%     fprintf('\n')
%     %
end
s0 = find(result(:,1)==min(result(:,1)));
Untitled2(sigma_rayvec(s0), result(s0,3), result(s0,2), n, [result(s0,4) result(s0,5)], method, result(s0,1))

%%
figure
plot(sigma_rayvec(1:10), result(:,1))
title('\sigma_{ray} Vs. Maximum (%) error')

figure
scatter(sigma_rayvec(1:10), result(:,2))
hold
coeff_rad_sqr = polyfit(sigma_rayvec(1:7)',result(1:7,2),1);
plot(x, coeff_rad_sqr(1).*x+coeff_rad_sqr(2),'r')
title('\sigma_{ray} Vs. radius (mm)')

figure
scatter(sigma_rayvec(1:10), result(:,3))
hold
coeff_sig_sqr = polyfit(sigma_rayvec(1:9)',result(1:9,3),1);
plot(x, coeff_sig_sqr(1).*x+coeff_sig_sqr(2),'r')
title('\sigma_{ray} Vs. \sigma_{sub}')

figure
plot(sigma_rayvec(1:10), result(:,6))
title('\sigma_{ray} Vs. computation time')

figure
scatter(sigma_rayvec(1:10), result(:,4))
hold
coeff= polyfit(sigma_rayvec(1:7)',result(1:7,4),2);
plot(x, coeff(1).*x.^2+coeff(2).*x+coeff(3),'r')
title('\sigma_{ray} Vs. normalization of weights')

figure
scatter(sigma_rayvec(1:10), result(:,5))
hold
coeff= polyfit(sigma_rayvec(1:9)',result(1:9,5),1);
plot(x, coeff(1).*x+coeff(2),'r')
title('\sigma_{ray} Vs. \sigma_w')
% 
coeff_rad_sqr3_correct = polyfit(sigma_rayvec(1:10)',result(:,2),1);
coeff_sig_sqr3_correct = polyfit(sigma_rayvec(1:10)',result(:,3),1);
coeff_sigW_sqr3_correct= polyfit(sigma_rayvec(1:9)',result(1:9,5),1);
coeff_w_sqr3_correct= polyfit(sigma_rayvec(1:7)',result(1:7,4),2);


%% Data analysis
% I decided not to use the last points of the data that I showed before.
% My decision comes form the fact that the error becomes bigger in those
% points (\sigma_{ray}>10) and this can be caused by the number of sub-rays
% that can be not sufficient. Further analysis will be done with n=3.
% Knowing this parameters we can predict the best sigma for the weighting
% gaussian and the best distance between the sub-ray shells.
% 
x = 3:0.01:12;
y = 8:0.01:17;

figure
hold
plot(x, coeff_sig_sqr_correct(1).*x+coeff_sig_sqr_correct(2),'r')
plot(y, coeff_sig_sqr3_correct(1).*y+coeff_sig_sqr3_correct(2),'b')
title('\sigma_{sub} square')
legend(strcat('n=2   ->   y=',num2str(coeff_sig_sqr_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sig_sqr_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_sig_sqr3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sig_sqr3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_sig_circ_correct(1).*x+coeff_sig_circ_correct(2),'r')
plot(y, coeff_sig_circ3_correct(1).*y+coeff_sig_circ3_correct(2),'b')
title('\sigma_{sub} circle')
legend(strcat('n=2   ->   y=',num2str(coeff_sig_circ_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sig_circ_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_sig_circ3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sig_circ3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_rad_sqr_correct(1).*x+coeff_rad_sqr_correct(2),'r')
plot(y, coeff_rad_sqr3_correct(1).*y+coeff_rad_sqr3_correct(2),'b')
title('r square')
legend(strcat('n=2   ->   y=',num2str(coeff_rad_sqr_correct(1),'%3.2f'),'x  +  ',num2str(coeff_rad_sqr_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_rad_sqr3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_rad_sqr3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_rad_circ_correct(1).*x+coeff_rad_circ_correct(2),'r')
plot(y, coeff_rad_circ3_correct(1).*y+coeff_rad_circ3_correct(2),'b')
title('r circle')
legend(strcat('n=2   ->   y=',num2str(coeff_rad_circ_correct(1),'%3.2f'),'x  +  ',num2str(coeff_rad_circ_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_rad_circ3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_rad_circ3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_sigW_sqr_correct(1).*x+coeff_sigW_sqr_correct(2),'r')
plot(y, coeff_sigW_sqr3_correct(1).*y+coeff_sigW_sqr3_correct(2),'b')
title('\sigma_w square')
legend(strcat('n=2   ->   y=',num2str(coeff_sigW_sqr_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sigW_sqr_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_sigW_sqr3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sigW_sqr3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_sigW_circ_correct(1).*x+coeff_sigW_circ_correct(2),'r')
plot(y, coeff_sigW_circ3_correct(1).*y+coeff_sigW_circ3_correct(2),'b')
title('\sigma_w circ')
legend(strcat('n=2   ->   y=',num2str(coeff_sigW_circ_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sigW_circ_correct(2),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_sigW_circ3_correct(1),'%3.2f'),'x  +  ',num2str(coeff_sigW_circ3_correct(2),'%3.2f')))

figure
hold
plot(x, coeff_w_sqr_correct(1).*x.^2+coeff_w_sqr_correct(2).*x+coeff_w_sqr_correct(3),'r')
plot(y, coeff_w_sqr3_correct(1).*y.^2+coeff_w_sqr3_correct(2).*y+coeff_w_sqr3_correct(3),'b')
title('normalization of weights square')
legend(strcat('n=2   ->   y=',num2str(coeff_w_sqr_correct(1),'%3.2f'),'x^2  +  ',num2str(coeff_w_sqr_correct(2),'%3.2f'), 'x  +  ',num2str(coeff_w_sqr_correct(3),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_w_sqr3_correct(1),'%3.2f'),'x^2  +  ',num2str(coeff_w_sqr3_correct(2),'%3.2f'), 'x  +  ',num2str(coeff_w_sqr3_correct(3),'%3.2f')))

figure
hold
plot(x, coeff_w_circ_correct(1).*x.^2+coeff_w_circ_correct(2).*x+coeff_w_circ_correct(3),'r')
plot(y, coeff_w_circ3_correct(1).*y.^2+coeff_w_circ3_correct(2).*y+coeff_w_circ3_correct(3),'b')
title('normalization of weights circle')
legend(strcat('n=2   ->   y=',num2str(coeff_w_circ_correct(1),'%3.2f'),'x^2  +  ',num2str(coeff_w_circ_correct(2),'%3.2f'), 'x  +  ',num2str(coeff_w_circ_correct(3),'%3.2f')),...
    strcat('n=3   ->   y=',num2str(coeff_w_circ3_correct(1),'%3.2f'),'x^2  +  ',num2str(coeff_w_circ3_correct(2),'%3.2f'), 'x  +  ',num2str(coeff_w_circ3_correct(3),'%3.2f')))



%% Results
% In general, we can say that the circular shape give us a better result
% (smaller maximum percentage error) than square shape. That is due to the
% circular simmetry of the problem. Square simmetry give us instead the
% possibility of a simplier construction of the grid of sub-beam and of
% first look maybe a way to save computaion time. More work will be done in
% this direction.
sigma_ray = 5;
n = 2;

[finalWeight, X1, sigma_sub, radius, posx, posy] = matRad_calcWeights(sigma_ray, n, 'square');
Untitled2(sigma_ray,sigma_sub,radius,n,X1,'square',1)

[finalWeight, X1, sigma_sub, radius, posx, posy] = matRad_calcWeights(sigma_ray, n, 'square');
Untitled2(sigma_ray,sigma_sub,radius,n,X1,'square',1)











%%
% % the comparison between the percentage error and the number of loops;
% % the comparison between the computation time and the number of loops;
% % load('weighted_data.mat');
% figure
% plot(pointi.*100,maxi)
% title('% Vs. #loops'); xlabel('#loops'); ylabel('% error');
% 
% figure
% plot(pointi.*100,timei)
% title('computation time Vs. #loops'); xlabel('#loops'); ylabel('seconds');
% %%
% % the comparison between the percentage error and the number of iterations;
% % the comparison between the computation time and the number of iterations;
% figure
% plot(pointi_iter.*10,maxi_iter)
% title('% Vs. #iter'); xlabel('#iter'); ylabel('% error');
% 
% figure
% plot(pointi_iter.*10,timei_iter)
% title('computation time Vs. #iter'); xlabel('#iter'); ylabel('seconds');
% 
% 
% 
% 
% 
% 

%%
% X1 = mR_findWeights(7, 6, 2, 3.05, 'circle');
% Untitled2(7,5.842,2.856,2,X1,'circle')
% 
% X1 = mR_findWeights(7, 5.684, 2, 2.856, 'circle');
% Untitled2(7,5.684,2.856,2,X1,'circle')
