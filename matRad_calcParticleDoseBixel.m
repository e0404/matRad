function bixel = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData, heteroCorrDepths, heteroCorrType, modulateBioDose, vTissueIndex)
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   dose = matRad_calcParticleDoseBixel(radDepths,radialDist_sq,SSD,focusIx,baseData)
%
% input
%   radDepths:      radiological depths
%   radialDist_sq:  squared radial distance in BEV from central ray
%   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
%   baseData:       base data required for particle dose calculation
%
% output
%   dose:   particle dose at specified locations as linear vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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

matRad_cfg =  MatRad_Config.instance();

% initialize heterogeneity correction
% function handle for calculating lateral dose
Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));

% function handle for calculating depth doses
sumGauss = @(x,mu,SqSigma,w) (1./sqrt(2*pi*ones(numel(x),1) .* SqSigma') .* ...
    exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) .* SqSigma' ))) * w;

% skip heterogeneity correction for other functions
if nargin < 5
    heteroCorrDepths = [];
end

% add potential offset
depths = baseData.depths + baseData.offset;

% convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
conversionFactor = 1.6021766208e-02;


%% interpolate depth dose, sigmas and weights and calculate lateral sigmas
if isstruct(baseData.Z) && ~isfield(baseData,'sigma1')
    
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,baseData.sigma,radDepths,'extrap');
    
    %compute lateral sigma
    sigmaSq = X.^2 + sigmaIni_sq;
    
elseif isstruct(baseData.Z) && isfield(baseData,'sigma1')
    
    % interpolate sigmas and weights
    X = matRad_interp1(depths,[baseData.sigma1 baseData.weight baseData.sigma2],radDepths);
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,1).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,3).^2 + sigmaIni_sq;
    
elseif ~isstruct(baseData.Z) && isfield(baseData,'sigma1')
    
    % interpolate depth dose, sigmas, and weights    
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma1 baseData.weight baseData.sigma2],radDepths,'extrap');
    
    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;
    
    % compute lateral sigmas
    sigmaSq_Narr = X(:,3).^2 + sigmaIni_sq;
    sigmaSq_Bro  = X(:,4).^2 + sigmaIni_sq;
  
else
    
    % interpolate depth dose and sigma
    X = matRad_interp1(depths,[conversionFactor*baseData.Z baseData.sigma],radDepths,'extrap');

    % set dose for query > tabulated depth dose values to zero
    X(radDepths > max(depths),1) = 0;
    
    %compute lateral sigma
    sigmaSq = X(:,2).^2 + sigmaIni_sq;
    
end

%  calculate lateral profiles
if isfield(baseData,'sigma1')
    
    L_Narr =  exp( -radialDist_sq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr); % Gauss
    L_Bro  =  exp( -radialDist_sq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
    
    bixel.L = baseData.LatCutOff.CompFac * ((1-X(:,2)).*L_Narr + X(:,2).*L_Bro); % (1-w)*L_Narr + w*L_Bro
    
else
    
    bixel.L = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) ./(2*pi*sigmaSq);
    
end

%% calculate sigma in range direction
if isstruct(baseData.Z)
    
    % calculate depthDoses with APM
    % no offset here...
    radDepths = radDepths - baseData.offset;
    
    % add sigma if heterogeneity correction wanted
    if ~isempty(heteroCorrDepths)
        switch heteroCorrType
            case 'complete'
                [~,lungDepthAtBraggPeakIx] = min(abs(radialDist_sq+(radDepths-baseData.peakPos).^2));
                lungDepthAtBraggPeak = heteroCorrDepths(lungDepthAtBraggPeakIx);
                ellSq = ones(numel(radDepths),1)* (baseData.Z.width'.^2 + matRad_getHeterogeneityCorrSigmaSq(lungDepthAtBraggPeak));
                
            case 'depthBased'
                % lungDepthAtGaussPeakIx = zeros(baseData.Z.mean)
                for i = 1:length(baseData.Z.mean)
                    [~,lungDepthAtGaussPeakIx(i)] = min(abs(radialDist_sq+(radDepths-baseData.Z.mean(i)).^2));
                end
                lungDepthAtGaussPeak = heteroCorrDepths(lungDepthAtGaussPeakIx);
                for i = 1:length(baseData.Z.mean)
                    ellSq(:,i) = ones(numel(radDepths),1)* (baseData.Z.width(i)'.^2 + matRad_getHeterogeneityCorrSigmaSq(lungDepthAtGaussPeak(i)));
                end
                
            case 'voxelwise'
                ellSq = bsxfun(@plus, baseData.Z.width'.^2, matRad_getHeterogeneityCorrSigmaSq(heteroCorrDepths));
                
            otherwise
                error('Error in heterogeneity correction')
        end
    else
        ellSq = ones(numel(radDepths),1)*baseData.Z.width'.^2;
    end
    
    % bixel.Z = (1./sqrt(2*pi*ellSq) .* exp(-bsxfun(@minus,baseData.Z.mean',radDepths).^2 ./ (2*ellSq)) )* baseData.Z.weight
    bixel.Z = sumGauss(radDepths,baseData.Z.mean,ellSq',baseData.Z.weight);
else
    
    bixel.Z = X(:,1);
    
    if ~isempty(heteroCorrDepths)
        warning('calcParticleDoseBixel: heterogeneity correction enabled but no APM base data was loaded.')
    end
    
end 
 
%% calculating the physical dose
bixel.physDose = conversionFactor * bixel.L .* bixel.Z;

%%
if ~isempty(heteroCorrDepths) && modulateBioDose
    if isfield(baseData,'alphaDose')
        % preallocate space for alpha beta
        tissueClasses = unique(vTissueIndex);
        bixel.Z_Aij = zeros(numel(radDepths),1);
        bixel.Z_Bij = zeros(numel(radDepths),1);
        
        for i = 1:numel(tissueClasses)
            
            ix = vTissueIndex == tissueClasses(i);
            bixel.Z_Aij(ix)  = conversionFactor * baseData.LatCutOff.CompFac * ...
                sumGauss(radDepths(ix)-baseData.offset,baseData.alphaDose(tissueClasses(i)).mean,...
                (baseData.alphaDose(tissueClasses(i)).width).^2,baseData.alphaDose(tissueClasses(i)).weight);
            
            bixel.Z_Bij(ix)  = conversionFactor * baseData.LatCutOff.CompFac * ...
                sumGauss(radDepths(ix)-baseData.offset,baseData.SqrtBetaDose(tissueClasses(i)).mean,...
                (baseData.SqrtBetaDose(tissueClasses(i)).width).^2,baseData.SqrtBetaDose(tissueClasses(i)).weight)';
            
        end
    else
        error('No alpha dose defined in base data, use fitted APM file with pln.propHeterogeneity.useDoseCurves')
    end
end

% check if we have valid dose values
if any(isnan(bixel.physDose)) || any(bixel.physDose<0)
   error('Error in particle dose calculation.');
end 
