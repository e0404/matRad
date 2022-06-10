function [robCube,robPassRate] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,slice,ct,cst,pln)
% robustness index calculation
%
% call
%   [robCube,robPassRate] = matRad_robustnessIndex(meanCube,stdCube,refDose,criteria,ct,cst,pln)
%
% input
%   meanCube:      dose cube as an M x N x O array
%   stdCube:       dose cube as an M x N x O array
%   criteria:      [1x2] vector specifying the dose robustness criteria;
%                  the first element is the acceptable percentage ratio between the difference
%                  of the calculated dose and the prescribed dose in TARGET with respect to the prescribed dose,
%                  second element is the acceptable percentage ratio between the std dose in TARGET and the prescribed dose
%   slice:         (optional) slice in cube1/2 that will be visualized
%   ct:            matRad ct structure
%   cst:           matRad cst structure
%   pln:           matRad pln structure
%
% output
%
%   robCube:       result of robustness index calculation
%   robPassRate:   rate of TARGET voxels passing the specified robustness criteria
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[env, ~] = matRad_getEnvironment();

% set parameters for robustness index calculation
if exist('criteria','var')
    meanDoseThreshold = criteria(1); % in [%]
    stdThreshold  = criteria(2); % in [%]
else
    meanDoseThreshold     = 5; % in [%]
    stdThreshold = 5; % in [%]
end

% set parameters for robustness index calculation
if ~exist('refDose','var') || refDose==0
    error('dose reference dose must be defined\n');
end

% check if cubes consistent
if ~isequal(size(meanCube),size(stdCube))
    error('dose cubes must be the same size\n');
end

if ~isfield(ct,'refScen')
    refScen=1;
else
    refScen=ct.refScen;
end

% Create target mask
targetMask = NaN(size(meanCube));
for  i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET')
        targetMask(cst{i,4}{refScen}) = 1;
    end
end

% Create dose mask
doseMask=(meanCube>0 | stdCube>0);

% evaluate robustness cube
meanDoseCrit = abs(meanCube-refDose)/refDose*100/meanDoseThreshold;
stdCrit  = stdCube/refDose*100/stdThreshold;
robCube = sqrt(meanDoseCrit.^2 + stdCrit.^2).*doseMask;

% compute robustness pass rate
numOfPassRobustness  = sum((robCube<=1).*~isnan(targetMask),'all');
robPassRate   = 100 * numOfPassRobustness / sum(~isnan(targetMask),'all');

% visualize if applicable
if exist('slice','var') && ~isempty(slice)
    
    plane      = 3;
    
    f=figure;
    f.Position(3:4) = [800 800];
    
    subplot(2,2,1);
    set(gcf,'Color',[1 1 1]);
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,stdCube/refDose,plane,slice,[],[],colorcube,[]);
    hold off;
    title('Relative uncertainty');
    
    subplot(2,2,2);
    set(gcf,'Color',[1 1 1]);
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,abs(meanCube-refDose*(meanCube>0))/refDose,plane,slice,[],[],colorcube,[]);
    hold off;
    title('Relative mean and prescription dose difference');
    
    maxRob = max(robCube,[],'all');
    doseWindow = [0 maxRob+0.1];
    
    mMap1=round(1/(maxRob)*256);
    mMap2=(256-mMap1);
    
    colormap1 = matRad_getColormap('priceOfRobustnessIndex',2*mMap1);
    colormap2 = matRad_getColormap('priceOfRobustnessIndex',2*mMap2);
    myColormap = [colormap1(1:mMap1,:); colormap2(mMap2+1:end,:)];
    
    subplot(2,2,3);
    set(gcf,'Color',[1 1 1]);
    
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,robCube,plane,slice,[],[],colorcube,myColormap,doseWindow);
    hold off;
    title('Robustness metric');

    subplot(2,2,4);
    set(gcf,'Color',[1 1 1]);
    
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,(robCube<=1).*doseMask,plane,slice,[],[],colorcube,[],[0 2.01]);
    hold off;
    title({[num2str(robPassRate,5) '% of points ' ...
        'pass robustness criterion (' num2str(meanDoseThreshold) '% / ' ...
        num2str(stdThreshold) '%)']});
    
    maxRob = max(robCube.*targetMask,[],'all');
    doseWindow = [0 maxRob+0.1];
    
    mMap1=round(1/(maxRob)*256);
    mMap2=(256-mMap1);
    
    colormap1 = matRad_getColormap('priceOfRobustnessIndex',2*mMap1);
    colormap2 = matRad_getColormap('priceOfRobustnessIndex',2*mMap2);
    myColormap = [colormap1(1:mMap1,:); colormap2(mMap2+1:end,:)];
    
    f = figure;
    numSlices = ct.cubeDim(3);
    %myColormap = matRad_getColormap('gammaIndex');
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,robCube.*targetMask,plane,slice,[],[],colorcube,myColormap,doseWindow);
    title({[num2str(robPassRate,5) '% of points ' ...
        'pass robustness criterion (' num2str(meanDoseThreshold) '% / ' ...
        num2str(stdThreshold) '%)']});
   
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback    = @(es,ed)  matRad_plotSliceWrapper2(gca,ct,cst,refScen,robCube.*targetMask,plane,round(es.Value),[],[],colorcube,myColormap,doseWindow);    
end
