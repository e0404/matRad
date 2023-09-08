matRad_rc;
matRad_cfg = MatRad_Config.instance();
matRad_cfg.propOpt.defaultMaxIter = 10000;
load 'TG119.mat'
%load('C:\r408i_data\r408i_data\CTDatasetMotion\102_HM10395_333.mat');

%load('PROSTATE.mat');
%cst{3,6} = [];
%% meta information for treatment plan (1) 
pln.numOfFractions  = 30;
pln.radiationMode   = 'protons';           % either photons / protons / helium / carbon
pln.machine         = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = [90 270]; % [?] ;
pln.propStf.couchAngles     = zeros(numel(pln.propStf.gantryAngles),1); % [?] ; 
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propDoseCalc.calcLET = 0;

pln.propOpt.runDAO          = false;      % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;      % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.propOpt.spatioTemp      = 0;
pln.propOpt.STscenarios     = 1;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

quantityOpt  = 'physicalDose';     % options: physicalDose, effect, RBExD
%=======================================> Model check error in bioModel
modelName    = 'none';             % none: for photons, protons, carbon            % constRBE: constant RBE for photons and protons 
                                   % MCN: McNamara-variable RBE model for protons  % WED: Wedenberg-variable RBE model for protons 
                                   % LEM: Local Effect Model for carbon ions


scenGenType  = 'wcScen';          % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen'                                          

% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation

pln.multScen = matRad_multScen(ct,scenGenType);
%pln.multScen.nSamples = 5;

%load('plnMultiScen.mat');
%% stf
stf = matRad_generateStf(ct,cst,pln);

%% cst

% for voiIdx=1:size(cst,1)
% 
%     for objIdx=1:size(cst{voiIdx,6},2)
%         if ~isempty(cst{voiIdx,6})
%             cst{voiIdx,6}{objIdx}.robustness = 'none'; 
%         end
%     end
% 
% 
% end
% method = 'STOCH';
% cst{1,6}{1}.robustness = method;
% cst{6,6}{1}.robustness = method;
% cst{7,6}{1}.robustness = method;
% cst{8,6}{1}.robustness = method;

% for voiIdx=1:size(cst,1)
%     if isequal(cst{voiIdx,3}, 'TARGET')
%         for objIdx=1:size(cst{voiIdx,6},2)
% 
%                 cst{voiIdx,6}{objIdx}.robustness = 'STOCH'; 
%         end
%     else
% 
%         for objIdx=1:size(cst{voiIdx,6},2)
% 
%                 cst{voiIdx,6}{objIdx}.robustness = 'none'; 
%         end
%     end
% end

for voiIdx=1:size(cst,1)
    if isequal(cst{voiIdx,3}, 'TARGET')
        %if ~isempty(cst{voiIdx,6})
            for objIdx=1:size(cst{voiIdx,6},2)
                cst{voiIdx,6}{objIdx}.robustness = 'STOCH';
            end
        %end
    else
        for objIdx=1:size(cst{voiIdx,6},2)

                cst{voiIdx,6}{objIdx}.robustness = 'none'; 
        end
    end
end

%% dij

pln.propDoseCalc.clearVoxelsForRobustness = 'none'; % none, targetOnly, oarOnly, objectivesOnly, [scenario indexes];


tic
dij_nominal  = matRad_calcParticleDose(ct,stf, pln,cst,0);
nominal_dij_time = toc;


tic
resultGUI_nominal = matRad_fluenceOptimization(dij_nominal,cst,pln);
nominal_opt_time = toc;




%% first method
%this is for excluding voxels a priori, brings numerical differences in the
%optimized weights
pln.propDoseCalc.clearVoxelsForRobustness = 'targetOnly'; % none, targetOnly, oarsOnly, objectivesOnly, [scenario indexes];

tic
%pln.propDoseCalc.clearMultiScenarioUnusedVoxels = true;

dij_reduced_first  = matRad_calcParticleDose(ct,stf, pln,cst,0);
reduced_dij_time_firstMethod = toc;


tic
resultGUI_reduced_first = matRad_fluenceOptimization(dij_reduced_first,cst,pln);
reduced_opt_time_firstMethod = toc;

w_diff_first = resultGUI_reduced_first.w - resultGUI_nominal.w;


%% second method test
%this computes the bixelDose regularly and only after sets voxels to zero.
%Brings no numerical difference 
pln.propDoseCalc.clearVoxelsForRobustness = 'targetOnly'; % none, targetOnly, oarsOnly, objectivesOnly, [scenario indexes];

tic
%pln.propDoseCalc.clearMultiScenarioUnusedVoxels = true;
dij_reduced  = matRad_calcParticleDose(ct,stf, pln,cst,0);
reduced_dij_time = toc;

tic
resultGUI_reduced = matRad_fluenceOptimization(dij_reduced,cst,pln);
reduced_opt_time = toc;
w_diff = resultGUI_reduced.w - resultGUI_nominal.w;

%% Gamma analysis
nominal_plan = dij_nominal.physicalDose{1}*resultGUI_nominal.w;
reduced_plan_first = dij_nominal.physicalDose{1}*resultGUI_reduced_first.w;
reduced_plan_second = dij_nominal.physicalDose{1}*resultGUI_reduced.w;

nominal_plan = reshape(nominal_plan, dij_nominal.doseGrid.dimensions);
reduced_plan_first = reshape(reduced_plan_first, dij_nominal.doseGrid.dimensions);
reduced_plan_second = reshape(reduced_plan_second, dij_nominal.doseGrid.dimensions);

slice = 31;
gamma_first = matRad_gammaIndex(nominal_plan,reduced_plan_first,[dij_nominal.doseGrid.resolution.x,dij_nominal.doseGrid.resolution.y,dij_nominal.doseGrid.resolution.z],[3 3]);
gamma_first = reshape(gamma_first, dij_nominal.doseGrid.dimensions);
cMap = matRad_getColormap('gammaIndex');


cst_n = matRad_resizeCstToGrid(cst,ct.x, ct.y, ct.z, dij_nominal.doseGrid.x, dij_nominal.doseGrid.y, dij_nominal.doseGrid.z);
f = figure;
imagesc(gamma_first(:,:,slice));
c.cubeDim = dij_nominal.doseGrid.dimensions;%{zeros(dij_nominal.doseGrid.dimensions)};
matRad_plotVoiContourSlice(gca(f), cst_n,c,1,1,3,slice);
colormap(cMap);
colorbar;
title('first method');

gamma_second= matRad_gammaIndex(nominal_plan,reduced_plan_second,[dij_nominal.doseGrid.resolution.x,dij_nominal.doseGrid.resolution.y,dij_nominal.doseGrid.resolution.z],[3 3]);
gamma_second = reshape(gamma_second, dij_nominal.doseGrid.dimensions);
cMap = matRad_getColormap('gammaIndex');


f = figure;
imagesc(gamma_second(:,:,slice));
c.cubeDim = dij_nominal.doseGrid.dimensions;%{zeros(dij_nominal.doseGrid.dimensions)};
matRad_plotVoiContourSlice(gca(f), cst_n,c,1,1,3,slice);
colormap(cMap);

colorbar;
title('second method');
%% Check
w = 1000*ones(size(dij_nominal.physicalDose{1},2),1);

nonEmptyScen = find(~cellfun(@isempty, dij_nominal.physicalDose));
targetIdx = cst_n{2,4}{1};

vi = [1:size(dij_nominal.physicalDose{1},1)]';
vi(targetIdx) = 0;

vi = vi(vi~=0);

for s=2:numel(nonEmptyScen)
     dist_nominal{s} = dij_nominal.physicalDose{nonEmptyScen(s)};%resultGUI_nominal.w;
     dist_reduced{s} = dij_reduced.physicalDose{nonEmptyScen(s)};%resultGUI_reduced.w;
     
     diffTotal{s} = dist_nominal{s}(targetIdx,:) - dist_reduced{s}(targetIdx,:);
     idxDiff{s} = find(diffTotal{s}>0);
     % diffTarget(s) = sum(dist_nominal{s}(targetIdx,:) - dist_reduced{s}(targetIdx,:), 'all');
     % diffOut(s) = sum(dist_nominal{s}(vi) - dist_reduced{s}(vi));
     % 
     % sumDistrib(s) = sum(dist_reduced{s});
     % sumDistribTarget(s) = sum(dist_reduced{s}(targetIdx));
     % 
     % zerosDist(s) = sum(dist_reduced{s}) - sum(dist_reduced{s}(targetIdx));
     % zeroDist2(s) = sumDistrib(s) -sumDistribTarget(s);
     % compl(s) = sum(dist_reduced{s}(vi));
end
%% asfsadf
 clear dist_nominal;
 clear dist_reduced;
 clear diffTarget;
 clear diffOut;
 clear sumDistrib;
 clear sumDistribTarget;
 clear zerosDist;
 clear zeroDist2;
 clear compl;

%% Comparison
for s=2:numel(nonEmptyScen)
     dist_nominal{s} = dij_nominal.physicalDose{nonEmptyScen(s)}*resultGUI_nominal.w;
     dist_reduced{s} = dij_reduced.physicalDose{nonEmptyScen(s)}*resultGUI_reduced.w;
end

 %% opt
pln.propOpt.clearUnusedVoxels = 0;

tic
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
originalTime = toc;

pln.propOpt.clearUnusedVoxels = 1;

tic
resultGUI_reduced = matRad_fluenceOptimization(dij,cst,pln);
afterTime = toc;

%% plots
figure;
subplot(1,2,1);
imagesc(resultGUI.physicalDose(:,:,80));
title('nominal plan');
subplot(1,2,2);
reducedDistribution = matRad_calcCubes(resultGUI_reduced.w, dij,1);
imagesc(reducedDistribution.physicalDose(:,:,80));
title('reduced plan');

diff = resultGUI.physicalDose - reducedDistribution.physicalDose;
max(max(max(diff)))
w_diff = resultGUI.w - resultGUI_reduced.w;
sum(w_diff)

%% Systematic

resolutions = [8,5,3,2];
originalTime = [];
originalIterations = [];
w_var = {};
w_var_tot = [];
afterTime = [];
afterIter = [];
for resolutionIdx = 1:size(resolutions,2)
    pln.propDoseCalc.doseGrid.resolution.x = resolutions(resolutionIdx);
    pln.propDoseCalc.doseGrid.resolution.y = resolutions(resolutionIdx);
    pln.propDoseCalc.doseGrid.resolution.z = resolutions(resolutionIdx);
    stf = matRad_generateStf(ct,cst,pln);
    dij  = matRad_calcParticleDose(ct,stf, pln,cst,0);

    pln.propOpt.clearUnusedVoxels = 0;

    tic
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    originalTime(resolutionIdx) = toc;
    originalIterations(resolutionIdx) = resultGUI.info.iter;

    pln.propOpt.clearUnusedVoxels = 1;
    
    tic
    resultGUI_reduced = matRad_fluenceOptimization(dij,cst,pln);
    afterTime(resolutionIdx) = toc;
    afterIter(resolutionIdx) = resultGUI_reduced.info.iter;

    w_var{resolutionIdx} = resultGUI.w - resultGUI_reduced.w;
    w_var_tot(resolutionIdx) = sum(w_var{resolutionIdx});
end




%% plots
figure;
plot(resolutions, originalTime, '.-');
hold on;
plot(resolutions, afterTime, '.-');
grid on;
grid minor;
legend('Original dij','Reduced dij');
xlabel('resolution [mm]', 'FontSize', 14);
ylabel('time[s]', 'FontSize', 14);
title('Prostate physicalDose', 'FontSize', 14);