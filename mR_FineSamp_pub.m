%% MatRad with Fine Sampling Beam
% In this document we will provide some examples on the results of MatRad
% implemented with Fine Sampling Beam Algorithm and some comparison with
% the previous version of MatRad.

%% Phantoms
% Here we show the phantoms on which we made our tests. The second one
% has a region of inhomogeneity in the center. 

addpath('tools')
load('BoxTest3_0_0.mat')
figure
subplot(1,2,1)
imagesc(ct2.cube{1}(:,:,80))
subplot(1,2,2)
imagesc(ct.cube{1}(:,:,80))

%% Results with homogeneous Phantom
% Shown after are the results for comparison between the algorithms in a
% homogeneous phantom for a single Ray

tic
matRad_calc_00 = matRad_calcDoseDirect(ct2,stf3,pln,cst,resultGUI.w);
t1 = toc
tic
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct2,stf3,pln,cst,resultGUI.w);
t2 = toc
figure
subplot(1,2,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
title(strcat('MatRad    time = ', num2str(t1)))
axis([65 95 30 100])
colorbar
subplot(1,2,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
axis([65 95 30 100])
colorbar
figure
subplot(1,2,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
title('Percentage difference')
axis([65 95 30 100])
colorbar
subplot(1,2,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
axis([65 95 30 100])
colorbar

%% Results with inhomogeneous Phantom
% Shown after are the results for comparison between the algorithms in an
% inhomogeneous phantom for a single Ray

tic
matRad_calc_00 = matRad_calcDoseDirect(ct,stf3,pln,cst,resultGUI.w);
t1 = toc
tic
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct,stf3,pln,cst,resultGUI.w);
t2 = toc
figure
subplot(1,2,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
title(strcat('MatRad    time = ', num2str(t1)))
axis([65 95 30 100])
colorbar
subplot(1,2,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
axis([65 95 30 100])
colorbar
figure
subplot(1,2,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
title('Percentage difference')
axis([65 95 30 100])
colorbar
subplot(1,2,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
axis([65 95 30 100])
colorbar

%% Results with inhomogeneous Phantom
% Shown after are the results for comparison between the algorithms in an
% inhomogeneous phantom for nine rays with square simmetry with central
% ray positioned on inhomogeneity position

matRad_calc_00 = matRad_calcDoseDirect(ct,stf2,pln,cst,resultGUI.w);
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct,stf2,pln,cst,resultGUI.w);
figure
subplot(1,2,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
title('MatRad')
axis([65 95 30 100])
colorbar
subplot(1,2,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title('MatRad Fine Sampling Beam')
axis([65 95 30 100])
colorbar
figure
subplot(1,2,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
title('Percentage difference')
axis([65 95 30 100])
colorbar
subplot(1,2,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
axis([65 95 30 100])
colorbar

%% Results with inhomogeneous Phantom at different angle
% Here we show the result for a beam (1 ray) entering the area with couch angle of
% 10 degrees and a gantry angle of 40 degrees in the previous seen phantoms.
% With and without inhomogeneities, respectevely.
% The red line represents the beginning of the inhomogeneity.

load('BoxTest3.mat')

tic
matRad_calc_00 = matRad_calcDoseDirect(ct2,stf6,pln,cst,resultGUI.w);
t1 = toc
tic
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct2,stf6,pln,cst,resultGUI.w);
t2 = toc
figure
subplot(2,1,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
title(strcat('MatRad    time = ', num2str(t1)))
axis([55 110 35 100])
colorbar
subplot(2,1,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
axis([55 110 35 100])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
title('Percentage difference')
axis([55 110 35 100])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
axis([55 110 35 100])
colorbar

%%
tic
matRad_calc_00 = matRad_calcDoseDirect(ct,stf4,pln,cst,resultGUI.w);
t1 = toc
tic
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct,stf4,pln,cst,resultGUI.w);
t2 = toc
figure
subplot(2,1,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
line([80 80],[20 120],'Color','r')
title(strcat('MatRad    time = ', num2str(t1)))
axis([55 110 45 100])
colorbar
subplot(2,1,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
line([80 80],[20 120],'Color','r')
title('Percentage difference')
axis([55 110 45 100])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar

%% Results with inhomogeneous Phantom at [90 0]
% Here we show the result for a beam (1 ray and 9 rays) entering the area with couch angle of
% 0 degrees and a gantry angle of 90 degrees in the previous seen phantom.
% The red line represents the beginning of the inhomogeneity.

load('BoxTest3_90_0.mat')

matRad_calc_00 = matRad_calcDoseDirect(ct,stf3,pln,cst,resultGUI.w);
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct,stf3,pln,cst,resultGUI.w);
figure
subplot(2,1,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
line([80 80],[20 120],'Color','r')
title('MatRad')
axis([55 110 45 100])
colorbar
subplot(2,1,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title('MatRad Fine Sampling Beam')
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
line([80 80],[20 120],'Color','r')
title('Percentage difference')
axis([55 110 45 100])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar

%%

matRad_calc_00 = matRad_calcDoseDirect(ct,stf2,pln,cst,resultGUI.w);
matRadFS_calc_00 = matRad_calcDoseDirectXXX(ct,stf2,pln,cst,resultGUI.w);
figure
subplot(2,1,1)
imagesc(matRad_calc_00.physicalDose(:,:,80))
line([80 80],[20 120],'Color','r')
title('MatRad')
axis([55 110 45 100])
colorbar
subplot(2,1,2)
imagesc(matRadFS_calc_00.physicalDose(:,:,80))
title('MatRad Fine Sampling Beam')
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_00.physicalDose(:,:,80)-matRadFS_calc_00.physicalDose(:,:,80)) ./max(max(matRad_calc_00.physicalDose(:,:,80))).*100)
line([80 80],[20 120],'Color','r')
title('Percentage difference')
axis([55 110 45 100])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_00.physicalDose,matRadFS_calc_00.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],80,1,'global',cst);
line([80 80],[20 120],'Color','r')
axis([55 110 45 100])
colorbar

%% Double wedges Phantom
% Here we have a figure representing the double wedges phantom on witch the
% experiment was run.

addpath('E:\Pezzano\MATLAB\Pezz\HIT Data\')
load('HITpoint_workspace.mat')
figure
imagesc(ct.cube{1}(:,:,round(pln.isoCenter(3)/10)))
text(50,250,'Water (1)','Color','k','FontSize',16)
text(200,150,'PMMA (1.165)','Color','k','FontSize',16)
text(400,400,'Air (0)','Color','w','FontSize',16)
arrow([450 250], [350 250], 'Beam')

%% Test on HIT double wedges Phantom for one sampling beam
% Here we compare the results for the simulation, made with MatRad, Syngo 
% and the new Fine Sampling Beam method, on the double wedges
% phantom for one sampling beam.

tic
matRad_calc_dw_p = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
t1 = toc
tic
matRadFS_calc_dw_p = matRad_calcDoseDirectXXX(ct,stf,pln,cst,resultGUI.w);
t2 = toc
figure
subplot(3,1,1)
imagesc(matRad_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title(strcat('MatRad    time = ', num2str(t1)))
axis([50 260 190 300])
colorbar
subplot(3,1,2)
imagesc(matRadFS_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
axis([50 260 190 300])
colorbar
subplot(3,1,3)
imagesc(resultGUI.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title('Syngo')
axis([50 260 190 300])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10))-resultGUI.physicalDose(:,:,round(pln.isoCenter(3)/10))) ./max(max(matRad_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10)))).*100)
title('Percentage difference between MatRad and Syngo')
axis([50 260 190 300])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_dw_p.physicalDose,resultGUI.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],round(pln.isoCenter(3)/10),1,'global',cst);
axis([50 260 190 300])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10))-matRadFS_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10))) ./max(max(matRad_calc_dw_p.physicalDose(:,:,round(pln.isoCenter(3)/10)))).*100)
title('Percentage difference between MatRad and FSB')
axis([50 260 190 300])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_dw_p.physicalDose,matRadFS_calc_dw_p.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],round(pln.isoCenter(3)/10),1,'global',cst);
axis([50 260 190 300])
colorbar


% save('HIT_dw_test1_point','matRad_calc_dw_p','matRadFS_calc_dw_p')

%% Test on HIT double wedges Phantom
% Here we compare the full simulation made with MatRad and the simulation made
% with the Fine Sampling Algorithm for the double wedges phantom.

load('HITcube_workspace.mat')
load('HIT_dw_test1.mat')

% tic
% matRad_calc_dw = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);
% t1 = toc
% tic
% matRadFS_calc_dw = matRad_calcDoseDirectXXX(ct,stf,pln,cst,resultGUI.w);
% t2 = toc
figure
subplot(3,1,1)
imagesc(matRad_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title(strcat('MatRad    time = ', num2str(t1)))
axis([50 260 190 300])
colorbar
subplot(3,1,2)
imagesc(matRadFS_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title(strcat('Fine Sampling Beam    time = ', num2str(t2)))
axis([50 260 190 300])
colorbar
subplot(3,1,3)
imagesc(resultGUI.physicalDose(:,:,round(pln.isoCenter(3)/10)))
title('Syngo')
axis([50 260 190 300])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10))-resultGUI.physicalDose(:,:,round(pln.isoCenter(3)/10))) ./max(max(matRad_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10)))).*100)
title('Percentage difference between MR and Syngo')
axis([50 260 190 300])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_dw.physicalDose,resultGUI.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],round(pln.isoCenter(3)/10),1,'global',cst);
axis([50 260 190 300])
colorbar
figure
subplot(2,1,1)
imagesc((matRad_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10))-matRadFS_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10))) ./max(max(matRad_calc_dw.physicalDose(:,:,round(pln.isoCenter(3)/10)))).*100)
title('Percentage difference between MR and FSB')
axis([50 260 190 300])
colorbar
subplot(2,1,2)
[gammaCube,gammaPassRateCell] = matRad_gammaIndex_NEW(matRad_calc_dw.physicalDose,matRadFS_calc_dw.physicalDose,[ct.resolution.x ct.resolution.y ct.resolution.z],[1 1],round(pln.isoCenter(3)/10),1,'global',cst);
axis([50 260 190 300])
colorbar

% save('HIT_dw_test1','matRad_calc_dw','matRadFS_calc_dw')


%%
% to be continued...








