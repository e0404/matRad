clear, clc

load protons_generic
load sqDev_generic.mat

% load protons_generic_TOPAS_cropped.mat
% load sqDev_topas.mat

spreadData = sqDev(4:17,:);
usedEnergies = [machine.data(1:6:end).energy];
usedEnergies = usedEnergies(4:17);
usedSpreads = linspace(0,3,102);


% spreadData = sqDev(1:15,:);
% usedEnergies = [machine.data(1:6:end).energy];
% usedSpreads = linspace(0,0.9,102);

sigEneSpr = zeros(3,numel(usedEnergies));

best = [];
for count = 1:size(usedEnergies,2)
    [~ ,i, ~] = intersect([machine.data(:).energy],usedEnergies(count));

%     p = polyfit(usedSpreads, spreadData(count,:), 2);
% 
%     best = [best, p(2)/(-2*p(1))];
%     
%     plot(usedSpreads, p(1)*usedSpreads.^2+p(2).*usedSpreads+p(3));
%     
%     hold on
%     scatter(usedSpreads, spreadData(count,:));
%     hold off
%     waitforbuttonpress;

    [~, spreadInd] = min(spreadData(count,:));
    sigEneSpr(3,count) = usedSpreads(spreadInd);
     
    newDepths = linspace(0,machine.data(i).depths(end),numel(machine.data(i).depths) * 100);
    newDose   = interp1(machine.data(i).depths, machine.data(i).Z, newDepths, 'spline');
     
    machine.data(i).depths = newDepths;
    machine.data(i).Z      = newDose;

    %interpolate range at 80% dose after peak.
    [maxV, maxI] = max(machine.data(i).Z);
    [~, r80ind] = min(abs(machine.data(i).Z(maxI:end) - 0.8 * maxV));
 

    %find FWHM w50 of bragg peak
    [~, d50rInd] = min(abs(machine.data(i).Z(maxI:end) - 0.5 * maxV));
    d50rInd = d50rInd - 1;
    d50_r = interp1(machine.data(i).Z(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                            machine.data(i).depths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
                        
    
         
    array1 = machine.data(i).depths(maxI-round(d50rInd):maxI+d50rInd);
    array2 = machine.data(i).Z(maxI-round(d50rInd):maxI+d50rInd);  
%     
%     plot(machine.data(i).depths, machine.data(i).Z);
%         d50rInd
% 
%     waitforbuttonpress;

%     array1 = machine.data(i).depths(maxI : end);
%     array2 = machine.data(i).Z(maxI :end);  


    gauss = @(x, a, sigma, b) a * exp(- (x - b).^2 / (2*sigma^2));

    funcs.objective = @(p) sum((gauss(array1, p(1), p(2),p(3)) - array2).^2);

    funcs.gradient = @(p) [ sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)));
                            sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* (array1 - p(3)).^2 / p(2)^3)
                            sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* 2 .* (array1 - p(3)) / (2 * p(2)^2))];


    options.lb = [0,  0, 0];
    options.ub = [ Inf, Inf, Inf];
    options.ipopt.limited_memory_update_type = 'bfgs';
            options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.print_level = 1;

    start = [maxV, 2 * d50_r, machine.data(i).depths(maxI)];
    [fitResult, ~] = ipopt (start, funcs, options);
    
%     plot(array1, gauss(array1, fitResult(1), fitResult(2),fitResult(3)));   
%     hold on
%     scatter(array1, array2)
%     hold off
%     fitResult(2);
% 
%     waitforbuttonpress;
    sigEneSpr(1,count) = fitResult(2);
    sigEneSpr(2,count) = machine.data(i).energy;
end
    
sigmas = sigEneSpr(1,:);
ene    = sigEneSpr(2,:);
best   = sigEneSpr(3,:);


