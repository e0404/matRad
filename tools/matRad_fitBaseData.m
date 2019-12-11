function fitData = matRad_fitBaseData(doseCube, resolution, energy, mcData, initSigma0, onAxis)
% fit analytical data to pencil beam stored in doseCube,
% pencil beam in positive y direction, target per default in center of
% plane
% 
% call
%    fitData = matRad_fitBaseData(doseCube, resolution, energy, initSigma0)
%
% input
%   doseCube:   dose cube as an M x N x O array
%   resolution: resolution of the cubes [mm/voxel]
%   energy:     energy of ray
%   initSigma0: initial sigma of beam, measured at entrance into doseCube
%   onAxis:     y and z coordinates of beam
%
% output 
%   fitData:    struct containing fit, structured in the same way as data
%               in machine.data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

% save cube dimesions
cubeDim = size(doseCube);


% set onAxis to default if not given
if ~exist('onAxis','var')
    onAxis = [cubeDim(2)/2 * resolution.y, cubeDim(3)/2 * resolution.z];
end

% extract IDD and discard values = 0
IDD = reshape(sum(sum(doseCube,2),3),cubeDim(1),1);
IDD = [IDD(1); IDD];
IDDnotZero = find(IDD);
IDD = IDD(IDDnotZero);

% calculate radiological depths for IDD
depthsIDD = [0, (resolution.x / 2 : resolution.x : resolution.x * cubeDim(1))];
depthsIDD = depthsIDD(IDDnotZero);

% interpolate IDD in steps of 0.05mm
IDD = interp1(depthsIDD, IDD, 0:0.05:depthsIDD(end), 'spline');
depthsIDD = 0:0.05:depthsIDD(end);

% calculate lateral radii to ray axis
axisGauss = zeros(cubeDim(2), cubeDim(3));
for i = 1:cubeDim(2)
    for j = 1:cubeDim(3)
       axisGauss(i,j) = sqrt((resolution.y * i - onAxis(1))^2 + (resolution.z * j - onAxis(2))^2); 
    end
end

resSigma1 = zeros(numel(IDDnotZero) - 1, 1);
resSigma2 = zeros(numel(IDDnotZero) - 1, 1);
resWeight = zeros(numel(IDDnotZero) - 1, 1);

preSigma2 = 20;

for i = 1:size(resSigma1,1)
% curve fitting gaussians

    % save y-z dose slice in profile
    profile = doseCube(i,:,:);
    profile = reshape(profile, cubeDim(2), cubeDim(3));
    profile = profile / sum(sum(profile)) / (resolution.x * resolution.y);    
    
    gauss1 = @(x , sigma1) 1 / (2 * pi * sigma1^2) * exp(- x.^2 / (2*sigma1^2));
    
    gauss2 = @(x, w, sigma1, sigma2) (1 - w) / (2 * pi * sigma1^2) .* exp(- x.^2 / (2*sigma1^2)) + ...
                                        w   / (2 * pi * sigma2^2) .* exp(- x.^2 / (2*sigma2^2)); 

    % first single gauss fit
    funcs1.objective = @(p) sum((gauss1(axisGauss, p(1)) - profile).^2, 'all');
    
    funcs1.gradient  = @(p) 2 * sum((gauss1(axisGauss, p(1)) - profile) .* ...  
                                exp(-axisGauss.^2 ./ (2 * p(1)^2)) .* (axisGauss.^2 - 2 .* p(1)^2) ./ (2 .* pi .* p(1)^5), 'all');
                             
                                     
    options1.lb =   0;
    options1.ub = Inf;
    options1.ipopt.hessian_approximation = 'limited-memory';
    options1.ipopt.limited_memory_update_type = 'bfgs';
    options1.ipopt.print_level = 1;
    start1 = 8;
    
    [fitResult1, ~] = ipopt (start1, funcs1, options1);                   
    preSigma1 = fitResult1;   
    
    %define fit parameters
    funcs2.objective = @(p) sum(sum((gauss2(axisGauss, p(1), p(2), p(3)) - profile).^2));

    funcs2.gradient = @(p) [ 2 * sum(sum((gauss2(axisGauss, p(1), p(2), p(3)) - profile) .* ...  
                                    (-1 / (2 * pi * p(2)^2) .* exp(-axisGauss.^2 / (2 * p(2)^2)) + 1 / (2 * pi * p(3)^2) .* exp(-axisGauss.^2 / (2 * p(3)^2)))));
                             2 * sum(sum((gauss2(axisGauss, p(1), p(2), p(3)) - profile) .* ...                                   
                                    (1 - p(1)) .* exp(-axisGauss.^2 / (2 * p(2)^2)) .* (axisGauss.^2 - 2 * p(2)^2) / (2 * pi * p(2)^5))); 
                             2 * sum(sum((gauss2(axisGauss, p(1), p(2), p(3)) - profile) .* ...
                                    p(1)       .* exp(-axisGauss.^2 / (2 * p(3)^2)) .* (axisGauss.^2 - 2 * p(3)^2) / (2 * pi * p(3)^5)))]; 

    options2.lb = [   0,   preSigma1 - 1,  preSigma1 + 1];
    options2.ub = [ 0.2,   preSigma1 + 1,  preSigma2 + 10];
%     options.ipopt.tol = 1e-30;

    options2.ipopt.hessian_approximation = 'limited-memory';
    options2.ipopt.limited_memory_update_type = 'bfgs';
    options2.ipopt.print_level = 1;

    %run fit and calculate actual sigma by squared substracting initial
    %sigma / spotsize

    start2 = [0.002, preSigma1, preSigma2];
    [fitResult, ~] = ipopt (start2, funcs2, options2);
    
    preSigma2 = fitResult(3);
    
    if (i == 1)
        if ~exist('initSigma0','var')
            noInitSigma = true;
            options.lb = 0;
            initSigma0 = 0;
        else
            noInitSigma = false;
            options.lb = initSigma0; 
        end
    end
    
    
    if (i == 1) && noInitSigma
        initSigma0 = fitResult(2);
    end
    
    
    if(fitResult(2) > initSigma0)
        actualSigma = sqrt(fitResult(2)^2 - initSigma0^2);
    else
        actualSigma = 0;
    end

    resSigma1(i) = actualSigma;
    resSigma2(i) = sqrt(fitResult(3)^2 - initSigma0^2);
    resWeight(i) = fitResult(1);
end   

% interpolate sigma on depths of IDD
resSigma1 = [0;            resSigma1];
resSigma2 = [resSigma2(1); resSigma2];
resWeight = [resWeight(1); resWeight];

depthsSigma = [0, (resolution.x / 2 : resolution.x : resolution.x * cubeDim(1))];
depthsSigma = depthsSigma(IDDnotZero);
resSigma1 = interp1(depthsSigma, resSigma1, depthsIDD);
resSigma2 = interp1(depthsSigma, resSigma2, depthsIDD);
resWeight = interp1(depthsSigma, resWeight, depthsIDD);
resSigma2 = smoothdata(resSigma2, 'gaussian', 200);

% interpolate range at 80% dose after peak.
[maxV, maxI] = max(IDD);
[~, r80ind] = min(abs(IDD(maxI:end) - 0.8 * maxV));
r80ind = r80ind - 1;
r80 = interp1(IDD(maxI + r80ind - 1:maxI + r80ind + 1), ...
                 depthsIDD(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

% conversion factor for IDDs
cf = 1 / 1.6021766208e-02 * resolution.y * resolution.z;

% save data in machine
fitData.energy =  energy;
fitData.range = r80;    
fitData.depths = depthsIDD';
fitData.Z = IDD' * cf;
fitData.peakPos = depthsIDD(maxI);
fitData.weight = resWeight';
fitData.offset = 0;

SAD = (2218 + 1839) / 2;
fitData.initFocus.dist  = linspace(SAD - 420, SAD + 420, 11);

corIso = (mcData.corNozzle * mcData.spotNozzle + mcData.divNozzle * mcData.z) / initSigma0;   
initSigmaZ = @(z) sqrt( initSigma0^2 - 2 * corIso * initSigma0 * mcData.divNozzle .* z + mcData.divNozzle^2 .* z.^2);


fitData.initFocus.sigma = initSigmaZ(-(fitData.initFocus.dist - SAD));
fitData.initFocus.SisFWHMAtIso = 2.3548 * initSigma0;

fitData.sigma1 = resSigma1';
fitData.sigma2 = resSigma2';


end
