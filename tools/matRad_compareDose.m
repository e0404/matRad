function [gammaCube,gammaPassRate,hfig] = matRad_compareDose(cube1, cube2, ct, cst, criteria, n,localglobal)
% Comparison of two dose cubes in terms of gamma index, absolute and visual difference
% 
% call
%    compareDose = matRad_compareDose(cube1,cube2,resolution,criteria,slice,interpoints,localglobal)
%
% input
%   cube1:         dose cube 1 as an M x N x O array
%   cube2:         dose cube 2 as an M x N x O array
%   ct:            plane view (coronal=1,sagittal=2,axial=3)
%   cst:           list of interesting volumes inside the patient
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
%%
resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];

%% Calculate iso-center slices
isoCenter=matRad_getIsoCenter(cst,ct,0);
slicename ={round(isoCenter(2)./resolution(2)),round(isoCenter(1)./resolution(1)),round(isoCenter(3)./resolution(3))};
doseWindowCube2 = [0 max([cube2(:); cube2(:)])];
doseWindowCube1 = [0 max([cube1(:); cube1(:)])];
planename={'coronal','sagittal','axial'};

%% Get the gamma cube
disp('Calculating gamma index cube...');
if exist('criteria','var')
    relDoseThreshold = criteria(1); % in [%]
    dist2AgreeMm     = criteria(2); % in [mm]
else
    dist2AgreeMm     = 3; % in [mm]
    relDoseThreshold = 3; % in [%]
end

[gammaCube,gammaPassRate] = matRad_gammaIndex(cube1,cube2,resolution,criteria,[],n,localglobal,cst);
set(gcf,'visible','off')


%%% Calculate absolute difference cube and dose windows for plots
DifferenceCube=abs(cube2-cube1);
doseWindowCube3 = [0 max([DifferenceCube(:); DifferenceCube(:)])];
gammadoseWindowCube = [0 max([gammaCube(:); gammaCube(:)])];

%% Plot everything
% Plot dose slices
for plane=1:3;
disp(['Plotting ',planename{plane},' plane...']);

% Initialize Figure
hfig.(planename{plane}).('fig')=figure('Renderer', 'painters', 'Position', [10 50 800 800]);
set(gcf,'Color',[1 1 1]);

% Plot Dose 1
hfig.(planename{plane}).('cube1').Axes=subplot(2,2,1);
    [hfig.(planename{plane}).('cube1').CMap,hfig.(planename{plane}).('cube1').Dose,hfig.(planename{plane}).('cube1').Ct,hfig.(planename{plane}).('cube1').Contour,hfig.(planename{plane}).('cube1').IsoDose]=matRad_plotSliceWrapper(gca,ct,cst,1,cube1,plane,slicename{plane},[],[],colorcube,[],doseWindowCube1,[],100);
    figtitle=get(gca,'title');
    figtitle=figtitle.String;

% Plot Dose 2
hfig.(planename{plane}).('cube2').Axes=subplot(2,2,2);
    [hfig.(planename{plane}).('cube2').CMap,hfig.(planename{plane}).('cube2').Dose,hfig.(planename{plane}).('cube2').Ct,hfig.(planename{plane}).('cube2').Contour,hfig.(planename{plane}).('cube2').IsoDose]=matRad_plotSliceWrapper(gca,ct,cst,1,cube2,plane,slicename{plane},[],[],colorcube,[],doseWindowCube2,[],100);

% Plot absolute difference
hfig.(planename{plane}).('diff').Axes=subplot(2,2,3);
    [hfig.(planename{plane}).('diff').CMap,hfig.(planename{plane}).('diff').Dose,hfig.(planename{plane}).('diff').Ct,hfig.(planename{plane}).('diff').Contour,hfig.(planename{plane}).('diff').IsoDose]=matRad_plotSliceWrapper(gca,ct,cst,1,DifferenceCube,plane,slicename{plane},[],[],colorcube,[],doseWindowCube3,[],100);
    colorbar;
    colormap(jet);

% Plot gamma analysis
hfig.(planename{plane}).('gamma').Axes=subplot(2,2,4);
    [hfig.(planename{plane}).('gamma').CMap,hfig.(planename{plane}).('gamma').Dose,hfig.(planename{plane}).('gamma').Ct,hfig.(planename{plane}).('gamma').Contour,hfig.(planename{plane}).('gamma').IsoDose]=matRad_plotSliceWrapper(gca,ct,cst,1,gammaCube,plane,slicename{plane},[],[],colorcube,[],gammadoseWindowCube,[],100);
    myColormap = matRad_getColormap('gammaIndex');
    colormap(gca,myColormap);
    colorbar

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

%% Visualize differences
% Through optimzation based on the biological effect we obtain a slightly
% different dose distribution as visualized by the following dose
% difference map

%if param.logLevel == 1
%    figure;
%    imagesc(resultGUI.RBExD(:,:,slice)-resultGUI_effect.RBExD(:,:,slice));
%    colorbar;
%    colormap(jet);
%end


%% 
% At this point we would like to see the absolute difference of the original optimization and the 
% recalculation. 

%if param.logLevel == 1
%    absDiffCube = resultGUI_effect.RBExD-resultGUI_tissue.RBExD;
%    figure,
%    matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);title('absolute difference')
%end

%% Plot gamma index

% matRad_gammaIndex();

disp('Done!');

end