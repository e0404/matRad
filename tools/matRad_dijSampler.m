function [cst,dij] = matRad_dijSampler(cst,ct,dij, marginSize,samplingFraction)
% sparce sampling of the dij for optimization
% 
% call
%   [cst,dij] = matRad_dijSampler(cst,ct,dij, marginSize,samplingFraction)
%
% input
%   dij:                matRad dij struct
%   cst:                matRad cst struct
%   ct:                 matRad ct struct
%   marginSize:         margin for each voi in pixel thickness
%   samplingFraction:   fraction of voxels to be sampled [ Target OAR]
%
% output
%   dij:                matRad dij struct
%   cst:                matRad cst struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = (2*marginSize)+1;
% s = strel('cube',m);
msk = ones(m,m,m);                  % mask for convolution

volThresh = 1000;                   % stuctures where the number of voxels is smaller that volThres will be ommitted from sampling
idx = [];
fullmap = [];                       % just for displaying
rng(1,'twister');                   % initialize random number generator (seed value )
cst = matRad_setOverlapPriorities(cst,ct.cubeDim);

for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        for j = 1:numel(cst{i,6})
            cst{i,6}(j).numOfVoxels = numel(cst{i,4}{1});
        end
    end
end

for i = 1:size(cst,1)
    
    % fraction of voxels sampled
    if (numel(cst{i,4}{1}) < volThresh)
        voiSamplingFraction = 1;
    elseif strcmp(cst{i,3}, 'TARGET')                   % threshold volume for small structures
        voiSamplingFraction = samplingFraction(1);
    elseif strcmp(cst{i,3}, 'OAR')
        voiSamplingFraction = samplingFraction(2);  % sampling fraction
    end
    % getting the margins
    mVOI = zeros(ct.cubeDim) ;
    mVOI(cst{i,4}{1}) = 1;
    %
    mVc = convn(mVOI,msk,'same');
    mVd = mVc - mVOI.*sum(msk(:));      % finding closed edge voxels
    idx = find(mVd<0);                  % idx
    % Finding a different way of extracting the margins
%     mVd = mVOI - imerode(mVOI,s);
%     idx = find(mVd<0);  
    % randomly sampling within the margins
    mVOI(idx) = 0;
    list = find(mVOI);
    if ~isempty(list)
        idx2 = list(randi([1 numel(list)],floor(voiSamplingFraction * numel(list)),1));
    else
        idx2 = [];
    end
    cst{i,4}{2} = cst{i,4}{1};          % saving the original idx ... required for ray casting !!!
    cst{i,4}{1} = [idx; idx2];
    
    idxRm =  setdiff(cst{i,4}{2},cst{i,4}{1});
    
    fullmap = [ fullmap; idx; idx2];    % for plotting the map of selected voxels
    
end
vec = zeros(numel(mVOI),1);
vec (fullmap)  = 1;
selectionDiag = spdiags(vec, 0 , dij.ctGrid.numOfVoxels, dij.ctGrid.numOfVoxels ) ;

tic,  dij.physicalDose{1} = ( dij.physicalDose{1}' * selectionDiag)';toc

if isfield(dij,'mAlphaDose')
    dij.mAlphaDose{1} = (dij.mAlphaDose{1}' * selectionDiag)';
    dij.mSqrtBetaDose{1} = (dij.mSqrtBetaDose{1}' * selectionDiag)';
end
xx = zeros(ct.cubeDim);
xx(fullmap) = 1;
figure, imagesc(xx(:,:,80));
end


