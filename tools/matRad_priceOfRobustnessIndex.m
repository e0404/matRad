function [resultGUISampRob] = matRad_priceOfRobustnessIndex(resultGUISampRob,doseCubeNom,doseCubeRob,ct,cst,refGy,refVol,refScen,priceWindow,calcType,slice)
% robustness index calculation
%
% call
%   [robCube,robPassRate] = matRad_robustnessIndex(meanCube,stdCube,p,criteria,ct,cst,pln)
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

matRad_cfg = MatRad_Config.instance();

[env, ~] = matRad_getEnvironment();

if ~exist('refScen','var') || isempty(refScen)
    refScen=1;
end

% set parameters for robustness index calculation
if ~exist('priceWindow','var') || isempty(priceWindow)
    minPrice = 10;
    maxPrice = -10;
else
    minPrice = priceWindow(1); %
    maxPrice  = priceWindow(2); %
end

% check if cubes consistent
if ~isequal(size(doseCubeNom),size(doseCubeRob))
    error('dose cubes must be the same size\n');
end

if ~exist('refVol', 'var') || isempty(refVol)
    refVol = [2 5 50 95 98];
end

if ~exist('refGy', 'var') || isempty(refGy)
    refGy = floor(linspace(0,max([doseCubeNom(:)]),6)*10)/10;
end

% Create target mask
OARMask = NaN(size(doseCubeNom));
for  i = 1:size(cst,1)
    if isequal(cst{i,3},'OAR')
        OARMask(cst{i,4}{refScen}) = 1;
    end
end

for  i = 1:size(cst,1)
    if isequal(cst{i,3},'TARGET')
        OARMask(cst{i,4}{refScen}) = NaN;
    end
end

indices     = find(OARMask);
numOfVoxels = numel(indices);
voiPrint    = sprintf('');

switch calcType
    case 'relative'
        priceCube=(doseCubeRob-doseCubeNom)./doseCubeNom;
    case 'absolute'
        priceCube=(doseCubeRob-doseCubeNom);
end

priceCube(doseCubeNom<=0) = 0;
priceCube(priceCube<minPrice) = minPrice;
priceCube(priceCube>maxPrice) = maxPrice;

qiNom = struct;
% get Dose, dose is sorted to simplify calculations
doseInVoiNom    = sort(doseCubeNom(indices));

if ~isempty(doseInVoiNom)
    % easy stats
    qiNom.mean = mean(doseInVoiNom);
    qiNom.std  = std(doseInVoiNom);
    qiNom.max  = doseInVoiNom(end);
    qiNom.min  = doseInVoiNom(1);
    
    DX = @(x) matRad_interp1(linspace(0,1,numOfVoxels),doseInVoiNom,(100-x)*0.01);
    VX = @(x) numel(doseInVoiNom(doseInVoiNom >= x)) / numOfVoxels;
    
    % create VX and DX struct fieldnames at runtime and fill
    for runDX = 1:numel(refVol)
        qiNom.(strcat('D_',num2str(refVol(runDX)))) = DX(refVol(runDX));
        voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
    end
    voiPrint = sprintf('%s\n%27s',voiPrint,' ');
    for runVX = 1:numel(refGy)
        sRefGy = num2str(refGy(runVX),3);
        qiNom.(['V_' strrep(sRefGy,'.','_') 'Gy']) = VX(refGy(runVX));
        voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '],voiPrint,VX(refGy(runVX))*100);
    end
    
end

voiPrint = sprintf('%s - Nominal opt. Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
    'OAR',qiNom.mean,qiNom.std,qiNom.max,qiNom.min,' ');

%We do it this way so the percentages in the string are not interpreted as format specifiers
matRad_cfg.dispInfo('%s\n',voiPrint);

qiRob = struct;
% get Dose, dose is sorted to simplify calculations
doseInVoiRob    = sort(doseCubeRob(indices));

if ~isempty(doseInVoiRob)
    % easy stats
    qiRob.mean = mean(doseInVoiRob);
    qiRob.std  = std(doseInVoiRob);
    qiRob.max  = doseInVoiRob(end);
    qiRob.min  = doseInVoiRob(1);
    
    DX = @(x) matRad_interp1(linspace(0,1,numOfVoxels),doseInVoiRob,(100-x)*0.01);
    VX = @(x) numel(doseInVoiRob(doseInVoiRob >= x)) / numOfVoxels;
    
    % create VX and DX struct fieldnames at runtime and fill
    for runDX = 1:numel(refVol)
        qiRob.(strcat('D_',num2str(refVol(runDX)))) = DX(refVol(runDX));
        voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
    end
    voiPrint = sprintf('%s\n%27s',voiPrint,' ');
    for runVX = 1:numel(refGy)
        sRefGy = num2str(refGy(runVX),3);
        qiRob.(['V_' strrep(sRefGy,'.','_') 'Gy']) = VX(refGy(runVX));
        voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '],voiPrint,VX(refGy(runVX))*100);
    end
    
end

voiPrint = sprintf('%s - Robust opt. Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
    'OAR',qiRob.mean,qiRob.std,qiRob.max,qiRob.min,' ');

%We do it this way so the percentages in the string are not interpreted as format specifiers
matRad_cfg.dispInfo('%s\n',voiPrint);

priceIndex=struct;
priceIndex.mean=qiRob.mean/qiNom.mean;

for runDX = 1:numel(refVol)
    if(qiNom.(strcat('D_',num2str(refVol(runDX))))>0)
        priceIndex.(strcat('D_',num2str(refVol(runDX)))) = qiRob.(strcat('D_',num2str(refVol(runDX))))/qiNom.(strcat('D_',num2str(refVol(runDX))));
    end
end

for runVX = 1:numel(refGy)
    sRefGy = num2str(refGy(runVX),3);
    if(qiNom.(['V_' strrep(sRefGy,'.','_') 'Gy'])>0)
        priceIndex.(['V_' strrep(sRefGy,'.','_') 'Gy']) = qiRob.(['V_' strrep(sRefGy,'.','_') 'Gy'])/qiNom.(['V_' strrep(sRefGy,'.','_') 'Gy']);
    end
end

resultGUISampRob.priceOfRobustnesAnalysis.priceCube=priceCube.*OARMask;
resultGUISampRob.priceOfRobustnesAnalysis.priceOfRobustness=priceIndex;

%voiPrint = sprintf('%s - Price of robustness = %5.2f +/- %5.2f \n%27s', ...
%    'OAR',priceIndex,(qiRob.std/abs(qiRob.std)+qiNom.std/abs(qiNom.std))*priceIndex,' ');

voiPrint = sprintf('%s - Price of robustness = %5.2f \n%27s', ...
    'OAR',priceIndex.mean,' ');

%We do it this way so the percentages in the string are not interpreted as format specifiers
matRad_cfg.dispInfo('%s\n',voiPrint);

% visualize if applicable
if exist('slice','var') && ~isempty(slice)
    
    plane      = 3;
    
    f = figure;
    numSlices = ct.cubeDim(3);
    doseWindow = [priceWindow(1)-0.1 priceWindow(2)+0.1];
    
    mMap1=round((1-maxPrice/(maxPrice-minPrice))*256);
    mMap2=(256-mMap1);
    
    colormap1 = matRad_getColormap('priceOfRobustnessIndex',2*mMap1);
    colormap2 = matRad_getColormap('priceOfRobustnessIndex',2*mMap2);
    myColormap = [colormap1(1:mMap1,:); colormap2(mMap2+1:end,:)];
    matRad_plotSliceWrapper2(gca,ct,cst,refScen,priceCube.*OARMask,plane,slice,[],[],colorcube,myColormap,doseWindow,[],[],'Price of robustness index');
    ax=gca;
    ax.Colorbar.TickLabels{1}=['<=' num2str(priceWindow(1))];
    ax.Colorbar.TickLabels{end}=['>=' num2str(priceWindow(2))];
    
    b = uicontrol('Parent',f,'Style','slider','Position',[50,5,420,23],...
        'value',slice, 'min',1, 'max',numSlices,'SliderStep', [1/(numSlices-1) , 1/(numSlices-1)]);
    b.Callback    = @(es,ed)  matRad_plotSliceWrapper2(gca,ct,cst,refScen,priceCube.*OARMask,plane,round(es.Value),[],[],colorcube,myColormap,doseWindow,[],[],'Price of robustness index');
    
end

end
