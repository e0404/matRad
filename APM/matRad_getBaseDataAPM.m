function [mMean, mWeight, mWidth, mLatSigma, mLUT, mDepth, vOffset, vCutOff, mCov] = matRad_getBaseDataAPM(ct,stf,pln,vTissueIndex,depthComp)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function converts the machine base data file into matrices and outputs
% a reduced number of APM components to increase the performance of evaluating
% the superposition of Gaussian components.
%
% call
%   [mMean, mWeight, mWidth, mLatSigma, LUT, mDepth, vOffset] = matRad_getReducesCompAPM(pln,depthComp)
%
% input
%   pln:              matRad plan meta information struct
%   depthComp:        string defining the depth component (e.g. 'Z' or 'alphaDose')
%
%
% output
%   mMean:         means of gaussian components for each inital beam energy
%   mWeight:       weight of gaussian components for each inital beam energy
%   mWidth:        squared sigma of gaussian components for each inital beam energy
%   mLatSigma:   squared sigma of beam widening in water per inital beam energy
%   mLUT:          look up table that see which component is active in
%                  which depth
%
%   mDepth:        depth values for each initial beam energy
%   vOffset:       offset to account for passive beam elements per inital
%                  beam energy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visBool = 0;

% function handle for calculating depth doses
sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
    exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

SGGauss  = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));

% prepare gauss base data
load([pln.radiationMode '_' pln.machine]);
numComp = length(machine.data(1).(depthComp)(1).weight);

% define wich tissue alphaDose should be used in case of bio opt
ixBio = unique(vTissueIndex);

if ixBio == 0
    ixBio = 1;
end

if numel(ixBio)>1
    matRad_dispToConsole('matRad_getReducedCompAPM: multiple tissues are not supported right now',[],'error');
end

if ~isfield(machine.data(1),'sigma')
    doubleGauss = true;
else
    doubleGauss = false;
end

%% compute SSDs
ctScen  = 1;
stf     = matRad_computeSSD(stf,ct,ctScen);

%% determine lateral cutoff
cutOffLevel          = 0.99;
visBoolLateralCutOff = 0;
vCutOff = [];
for i = 1:size(stf,2)
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),ctScen,visBoolLateralCutOff);
    for j = 1:numel(machine.data)
        vCutOff{j,i} = machine.data(j).LatCutOff;
    end
    
end

% loop over all entries
for i = 1:length(machine.data)
    
    %square the sigmas
    machine.data(i).(depthComp)(ixBio).width = machine.data(i).(depthComp)(ixBio).width .^2;
    
    %% find reduced components
    if pln.probOpt.useReducedComp
        
        numDepths                       = length(machine.data(i).depths);
        mContributions                  = zeros(numComp,numDepths);
        machine.data(i).(depthComp)(ixBio).LUT = false(numComp,numDepths);
        
        if visBool; figure,hold on;  end
        % loop over the components
        
        for iComp = 1:numComp
            vY = machine.data(i).(depthComp)(ixBio).weight(iComp) * SGGauss(machine.data(i).depths,machine.data(i).(depthComp)(ixBio).mean(iComp),machine.data(i).(depthComp)(ixBio).width(iComp));
            mContributions(iComp,:) = vY;
            machine.data(i).(depthComp)(ixBio).LUT(iComp,:) = vY > (max(vY)*pln.probOpt.reducedComprelDose );  % include 99 percent of each gauss
            if visBool; plot(machine.data(i).depths,vY);end
        end
        
        machine.data(i).(depthComp)(ixBio).numDepth = numDepths;
        doseRef                              = sumGauss(machine.data(i).depths,machine.data(i).(depthComp)(ixBio).mean,machine.data(i).(depthComp)(ixBio).width,machine.data(i).(depthComp)(ixBio).weight);
        doseReduced                          = zeros(numel(machine.data(i).depths),1);
        
        for j = 1:numel(machine.data(i).depths)
            doseReduced(j)    = sumGauss(machine.data(i).depths(j),machine.data(i).(depthComp)(ixBio).mean((machine.data(i).(depthComp)(ixBio).LUT(:,j))),...
                machine.data(i).(depthComp)(ixBio).width(machine.data(i).(depthComp)(ixBio).LUT(:,j)),...
                machine.data(i).(depthComp)(ixBio).weight(machine.data(i).(depthComp)(ixBio).LUT(:,j)));
        end
        
        if visBool
            figure,hold on,
            plot(machine.data(i).depths,doseRef,'b','LineWidth',2);   plot(machine.data(i).depths,doseReduced,'r:','LineWidth',2), grid on, grid minor,
            legend({'org profile','reduced components'})
        end
        
    else
        machine.data(i).(depthComp)(ixBio).LUT = [];
    end
    
end


mMean     = zeros(numComp,numel(machine.data));                                           % mean of gaussian components
mWeight   = zeros(numComp,numel(machine.data));                                           % weight of gaussian components
mWidth    = zeros(numComp,numel(machine.data));                                           % squared sigma of gaussian components

mDepth    = NaN*(zeros(numel(machine.data(end).depths),numel(machine.data)));               % depth values (query points) of the depth dose profiles
mLUT      = logical((zeros(numel(machine.data),numel(machine.data(end).depths),numComp)));     % look up table to see which components are are 'activated' and deactivated
vOffset   = zeros(numel(machine.data),1);                                                   % offset vector to account for the rad depth. of passive beam line elements

mCov      = cell(length(machine.data),1);

% create matrix for lateral beam spread in water
if doubleGauss
    mLatSigma = NaN*(zeros(numel(machine.data(end).depths),numel(machine.data),3));
else
    mLatSigma = NaN*(zeros(numel(machine.data(end).depths),numel(machine.data)));
end

% loop over all entries
for i = 1:numel(machine.data)
    
    mMean(:,i)   = machine.data(i).(depthComp)(ixBio).mean;
    mWeight(:,i) = machine.data(i).(depthComp)(ixBio).weight;
    mWidth(:,i)  = machine.data(i).(depthComp)(ixBio).width;                                       % with is already squared
    mDepth(1:numel(machine.data(i).depths),i) = machine.data(i).depths;
    
    if doubleGauss
        mLatSigma(1:numel(machine.data(i).depths),i,1) = machine.data(i).sigma1;
        mLatSigma(1:numel(machine.data(i).depths),i,2) = machine.data(i).weight;
        mLatSigma(1:numel(machine.data(i).depths),i,3) = machine.data(i).sigma2;
        
    else
        mLatSigma(1:numel(machine.data(i).depths),i) = machine.data(i).sigma;
    end
    

    % save look up table
    if pln.probOpt.useReducedComp
        mLUT(i,1:numel(machine.data(i).depths),:) = logical(machine.data(i).(depthComp)(ixBio).LUT)';
    else
        mLUT(i,1:numel(machine.data(i).depths),:) = true;
    end
    
    vOffset(i) = machine.data(i).offset;
 
end

mLUT   = permute(mLUT,[3 2 1]);

mDepth = bsxfun(@plus,mDepth,vOffset');


end % end of function



