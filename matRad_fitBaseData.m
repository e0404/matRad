function fitData = matRad_fitBaseData(doseCube, resolution, energy, initSigma0, onAxis)
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

resSigma = zeros(numel(IDDnotZero) - 1, 1);
for i = 1:size(resSigma,1)
% curve fitting gaussians

    % save y-z dose slice in profile
    profile = doseCube(i,:,:);
    profile = reshape(profile, cubeDim(2), cubeDim(3));
    
    gauss = @(r, a, sigma) a * exp(- r.^2 / (2*sigma^2));

    % define fit parameters
    funcs.objective = @(p) sum(sum((gauss(axisGauss, p(1), p(2)) - profile).^2));

    funcs.gradient = @(p) [ sum(sum(2 * (p(1) * exp(-axisGauss.^2 / (2 * p(2)^2)) - profile) ...
                                            .* exp(-axisGauss.^2 / (2 * p(2)^2))));
                            sum(sum(2 * (p(1) * exp(-axisGauss.^2 / (2 * p(2)^2)) - profile) ...
                                            .* p(1) .* exp(-axisGauss.^2 / (2 * p(2)^2)) .* axisGauss.^2 / p(2)^3))];

    options.lb = [0,  0];
    options.ub = [ Inf,  Inf];

    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.limited_memory_update_type = 'bfgs';
    options.ipopt.print_level = 1;
    options.ipopt.tol = 1e-16;

    % run fit and calculate actual sigma by squared substracting initial
    % sigma
    start = [max(max(profile)); initSigma0];
    [fitResult, ~] = ipopt (start, funcs, options);

    if(fitResult(2) > initSigma0)
        actualSigma = sqrt(fitResult(2)^2 - initSigma0^2);
    else
        actualSigma = 0;
    end

    resSigma(i) = actualSigma;
end   

% interpolate sigma on depths of IDD
resSigma = [0; resSigma];
depthsSigma = [0, (resolution.x / 2 : resolution.x : resolution.x * cubeDim(1))];
depthsSigma = depthsSigma(IDDnotZero);
resSigma = interp1(depthsSigma, resSigma, depthsIDD);

% interpolate range at 80% dose after peak.
[maxV, maxI] = max(IDD);
[~, r80ind] = min(abs(IDD(maxI:end) - 0.8 * maxV));
r80ind = r80ind - 1;
r80 = interp1(IDD(maxI + r80ind - 1:maxI + r80ind + 1), ...
                 depthsIDD(maxI + r80ind - 1:maxI + r80ind + 1), 0.8 * maxV);

% conversion factor for IDDs
cf = 1 / 1.6021766208e-02 * resolution.y * resolution.z;

% save data in machine
fitData.range = r80;    
fitData.energy =  energy;
fitData.depths = depthsIDD';
fitData.Z = IDD' * cf;
fitData.peakPos = depthsIDD(maxI);
fitData.sigma = resSigma';
fitData.offset = 0;

fitData.initFocus.dist  = [0, 20000];
fitData.initFocus.sigma = [initSigma0, initSigma0];
fitData.initFocus.SisFWHMAtIso = 2.3548 * initSigma0;
end
