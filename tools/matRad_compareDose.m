function [gammaCube,gammaPassRate,hfig] = matRad_compareDose(cube1, cube2, ct, cst,enable , contours, pln, criteria, n,localglobal)
% Comparison of two dose cubes in terms of gamma index, absolute and visual difference
%
% call
%    compareDose = matRad_compareDose(cube1,cube2,ct,cst,criteria,n,localglobal)
%
% input
%   cube1:         dose cube 1 as an M x N x O array
%   cube2:         dose cube 2 as an M x N x O array
%   ct:            ct struct with ct cube
%   cst:           list of interesting volumes inside the patient as matRad
%                  struct (optional, does not calculate gamma and DVH)
%   enable         (optional) specify if all sections are evaluated
%                  boolean 3x1 array
%                  1:
%                  2:
%                  3:
%   contours       (optional) specify if contours are plotted, boolean
%   pln            (optional) specify BioModel for DVH plot
%   criteria:      [1x2] vector (optional) specifying the distance to agreement
%                  criterion; first element is percentage difference,
%                  second element is distance [mm], default [3 3]
%   n:             number of interpolations (optional). there will be 2^n-1
%                  interpolation points. The maximum suggested value is 3.
%                  default n=0
%   localglobal:   parameter to choose between 'global' and 'local'
%                  normalization (optional)
%
%
% output
%
%   gammaCube:      result of gamma index calculation
%   gammaPassRate:  rate of voxels passing the specified gamma criterion
%                   evaluated for every structure listed in 'cst'.
%                   note that only voxels exceeding the dose threshold are
%                   considered.
%   hfig:           Figure handle struct for all 3 figures and subplots,
%                   indexed by plane names
%
% References gamma analysis:
%   [1]  http://www.ncbi.nlm.nih.gov/pubmed/9608475
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

%% check if cubes consistent
if ~isequal(size(cube1),size(cube2))
    error('dose cubes must be the same size\n');
end

if ~exist('localglobal','var')
    localglobal = 'global';
end
if ~exist('n','var')
    n = 0;
end
if ~exist('criteria','var')
    criteria = [3 3];
end
if ~exist('cst','var') || isempty(cst)
    cst = [];
    fprintf('[\bWarning: cst not specified, skipping DVH, isocenter in the center of the doseCube]\b\n');
end
if ~exist('pln','var') || isempty(pln)
    pln = [];
end
if ~exist('contours','var') || isempty(contours)
    contours = true;
end

if (~exist('enable','var') || isempty(enable)) || strcmp(enable,'all')
    enable = [1 1 1];
end
if enable(1)==0;
    gammaCube=[];
    gammaPassRate=[];
end
%% Load colormap for difference
load('diffCMap.mat');

%% Calculate iso-center slices and resolution
if isempty(cst)
    s=size(cube1);
    isoCenter = [s(1)/2 s(2)/2 s(3)/2];
else
    isoCenter = matRad_getIsoCenter(cst,ct,0);
end

resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];

slicename = {round(isoCenter(2)./resolution(2)),round(isoCenter(1)./resolution(1)),round(isoCenter(3)./resolution(3))};
doseWindow = [0 max([cube1(:); cube2(:)])];
planename = {'coronal','sagittal','axial'};

%% Get the gamma cube
if enable(1)==1;
    disp('Calculating gamma index cube...');
    if exist('criteria','var')
        relDoseThreshold = criteria(1); % in [%]
        dist2AgreeMm     = criteria(2); % in [mm]
    else
        dist2AgreeMm     = 3; % in [mm]
        relDoseThreshold = 3; % in [%]
    end
    
    [gammaCube,gammaPassRate] = matRad_gammaIndex(cube1,cube2,resolution,criteria,[],n,localglobal,cst);
    
    
    %%% Calculate absolute difference cube and dose windows for plots
    differenceCube  = cube1-cube2;
    doseDiffWindow  = [-max(differenceCube(:)) max(differenceCube(:))];
    doseGammaWindow = [0 max(gammaCube(:))];
    
    
    %% Plot everything
    % Plot dose slices
    if contours == false
        cst = [];
    end
    
    for plane=1:3;
        disp(['Plotting ',planename{plane},' plane...']);
        
        % Initialize Figure
        hfig.(planename{plane}).('fig') = figure('Renderer', 'painters', 'Position', [10 50 800 800]);
        set(gcf,'Color',[1 1 1]);
        
        % Plot Dose 1
        hfig.(planename{plane}).('cube1').Axes = subplot(2,2,1);
        [hfig.(planename{plane}).('cube1').CMap,...
            hfig.(planename{plane}).('cube1').Dose,...
            hfig.(planename{plane}).('cube1').Ct,...
            hfig.(planename{plane}).('cube1').Contour,...
            hfig.(planename{plane}).('cube1').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slicename{plane},[],[],colorcube,jet,doseWindow,[],100);
        figtitle = get(gca,'title');
        figtitle = figtitle.String;
        
        % Plot Dose 2
        hfig.(planename{plane}).('cube2').Axes = subplot(2,2,2);
        [hfig.(planename{plane}).('cube2').CMap,...
            hfig.(planename{plane}).('cube2').Dose,...
            hfig.(planename{plane}).('cube2').Ct,...
            hfig.(planename{plane}).('cube2').Contour,...
            hfig.(planename{plane}).('cube2').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slicename{plane},[],[],colorcube,jet,doseWindow,[],100);
        
        % Plot absolute difference
        hfig.(planename{plane}).('diff').Axes = subplot(2,2,3);
        [hfig.(planename{plane}).('diff').CMap,...
            hfig.(planename{plane}).('diff').Dose,...
            hfig.(planename{plane}).('diff').Ct,...
            hfig.(planename{plane}).('diff').Contour,...
            hfig.(planename{plane}).('diff').IsoDose] = ...
            matRad_plotSliceWrapper(gca,ct,cst,1,differenceCube,plane,slicename{plane},[],[],colorcube,diffCMap,doseDiffWindow,[],100);
        
        % Plot gamma analysis
        hfig.(planename{plane}).('gamma').Axes = subplot(2,2,4);
        gammaCMap = matRad_getColormap('gammaIndex');
        [hfig.(planename{plane}).('gamma').CMap,...
            hfig.(planename{plane}).('gamma').Dose,...
            hfig.(planename{plane}).('gamma').Ct,...
            hfig.(planename{plane}).('gamma').Contour,...
            hfig.(planename{plane}).('gamma').IsoDose]=...
            matRad_plotSliceWrapper(gca,ct,cst,1,gammaCube,plane,slicename{plane},[],[],colorcube,gammaCMap,doseGammaWindow,[],100);
        
        set(hfig.(planename{plane}).('fig'),'name',figtitle);
        
        %% Adjusting axes
        
        matRad_plotAxisLabels(hfig.(planename{plane}).('cube1').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube1').Axes, 'title'), 'string', 'Dose 1');
        matRad_plotAxisLabels(hfig.(planename{plane}).('cube2').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('cube2').Axes, 'title'), 'string', 'Dose 2');
        matRad_plotAxisLabels(hfig.(planename{plane}).('diff').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('diff').Axes, 'title'), 'string', 'Absolute difference');
        matRad_plotAxisLabels(hfig.(planename{plane}).('gamma').Axes,ct,plane,slicename{plane},[],100);
        set(get(hfig.(planename{plane}).('gamma').Axes, 'title'), 'string', {[num2str(gammaPassRate{1,2},5) '% of points > ' num2str(relDoseThreshold) '% pass gamma criterion (' num2str(relDoseThreshold) '% / ' num2str(dist2AgreeMm) 'mm)']; ['with ' num2str(2^n-1) ' interpolation points']});
        
    end
end
%% Plot profiles through isoCenter

if enable(2)==1;
    disp(['Plotting profiles...']);
    fontsize=12;
    profilex{1} = cube1(:,slicename{2},slicename{3});
    profiley{1} = permute(cube1(slicename{1},slicename{2},:),[3 2 1]);
    profilez{1} = permute(cube1(slicename{1},:,slicename{3}),[2 3 1]);
    
    profilex{2} = cube2(:,slicename{2},slicename{3});
    profiley{2} = permute(cube2(slicename{1},slicename{2},:),[3 2 1]);
    profilez{2} = permute(cube2(slicename{1},:,slicename{3}),[2 3 1]);
    
    cubes={'cube1','cube2','Dose 1','Dose 2'};
    %%
    if exist('pln','var') && ~isempty(pln)
        if strcmp(pln.bioParam.model,'none')
            l='Dose [Gy]';
        else
            l='RBE x Dose [Gy(RBE)]';
        end
    else
        l='Dose [Gy]';
    end
    
    %%
    for i=1:2;
        currentCube=cubes{i};
        hfig.profile.(currentCube).fig = figure('Renderer', 'painters', 'Position', [10 50 800 800]);
        set(gcf,'Color',[1 1 1]);
        
        hfig.profile.(currentCube).x = subplot(2,2,1);
        plot(profilex{i})
        xlabel('X [mm]','FontSize',fontsize)
        ylabel(l,'FontSize',fontsize);
        title(['x-Profile ',cubes{i+2}]);
        
        hfig.profile.(currentCube).y = subplot(2,2,2);
        plot(profiley{i})
        xlabel('Y [mm]','FontSize',fontsize)
        ylabel(l,'FontSize',fontsize);
        title(['y-Profile ',cubes{i+2}]);
        
        hfig.profile.(currentCube).z = subplot(2,2,3);
        plot(profilez{i})
        xlabel('Z [mm]','FontSize',fontsize)
        ylabel(l,'FontSize',fontsize);
        title(['z-Profile ',cubes{i+2}]);
        
        set(hfig.profile.(currentCube).fig,'name',['Profiles ',cubes{i+2},':, x=',num2str(slicename{1}),'mm, y=',num2str(slicename{2}),'mm, z=',num2str(slicename{3}),'mm']);
    end
end
%% Calculate and plot DVH

if enable(3)==1;
    if ~isempty(cst)
        disp('Calculating DVH...');
        dvh1=matRad_calcDVH(cst,cube1);
        dvh2=matRad_calcDVH(cst,cube2);
        % Plot DVH
        disp(['Plotting DVH...']);
        
        hfig.dvh.cube1 = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
        set(gcf,'Color',[1 1 1]);
        matRad_showDVH(dvh1,cst,pln)
        title('Dose 1')
        
        hfig.dvh.cube2 = figure('Renderer', 'painters', 'Position', [10 100 1000 700]);
        set(gcf,'Color',[1 1 1]);
        matRad_showDVH(dvh2,cst,pln)
        title('Dose 2')
    end
end
%%
disp('Done!');

end