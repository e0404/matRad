%% Example: Generate your own phantom geometry
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% In this example we will show how to calculate the variance of a carbon
% ion treatment plan and how to perform probabilistic optimization by
% optimization the expectation value of a squared deviation objective
clc, clear, close all

%% Create a CT image series
xDim = 130;
yDim = 130;
zDim = 50;

ct.cubeDim      = [xDim yDim zDim];
ct.resolution.x = 1.7;
ct.resolution.y = 1.7;
ct.resolution.z = 2;
ct.numOfCtScen  = 1;
 
% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target

ixOAR = 1;
ixPTV = 2;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';
 
% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 2;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,6}.type        = 'square overdosing';
cst{ixOAR,6}.dose        = 20;
cst{ixOAR,6}.penalty     = 50;
cst{ixOAR,6}.EUD         = NaN;
cst{ixOAR,6}.volume      = NaN;
cst{ixOAR,6}.coverage    = NaN;
cst{ixOAR,6}.robustness  = 'none';

cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,6}.type        = 'square deviation';
cst{ixPTV,6}.dose        = 60;
cst{ixPTV,6}.penalty     = 300;
cst{ixPTV,6}.EUD         = NaN;
cst{ixPTV,6}.volume      = NaN;
cst{ixPTV,6}.coverage    = NaN;
cst{ixPTV,6}.robustness  = 'none';


%% Lets create either a cubic or a spheric phantom
TYPE = 'spheric';   % either 'cubic' or 'spheric'

% first the OAR
cubeHelper = zeros(ct.cubeDim);

switch TYPE
   
   case {'cubic'}
      
      xLowOAR  = round(xDim/2 - xDim/4);
      xHighOAR = round(xDim/2 + xDim/4);
      yLowOAR  = round(yDim/2 - yDim/4);
      yHighOAR = round(yDim/2 + yDim/4);
      zLowOAR  = round(zDim/2 - zDim/4);
      zHighOAR = round(zDim/2 + zDim/4);
      
      for x = xLowOAR:1:xHighOAR
         for y = yLowOAR:1:yHighOAR
            for z = zLowOAR:1:zHighOAR
               cubeHelper(x,y,z) = 1;
            end
         end
      end
      
   case {'spheric'}
      
      radiusOAR = xDim/6;
      
      for x = 1:xDim
         for y = 1:yDim
            for z = 1:zDim
               currPost = [x y z] - round([ct.cubeDim./2]);
               if  sqrt(sum(currPost.^2)) < radiusOAR
                  cubeHelper(x,y,z) = 1;
               end
            end
         end
      end
      
end

% extract the voxel indices and save it in the cst
cst{ixOAR,4}{1} = find(cubeHelper);


% second the PTV
cubeHelper = zeros(ct.cubeDim);

switch TYPE
   
   case {'cubic'}
      
      xLowPTV  = round(xDim/2 - xDim/8);
      xHighPTV = round(xDim/2 + xDim/8);
      yLowPTV  = round(yDim/2 - yDim/8);
      yHighPTV = round(yDim/2 + yDim/8);
      zLowPTV  = round(zDim/2 - zDim/8);
      zHighPTV = round(zDim/2 + zDim/8);
      
      cubeHelper = zeros(ct.cubeDim);
      
      for x = xLowPTV:1:xHighPTV
         for y = yLowPTV:1:yHighPTV
            for z = zLowPTV:1:zHighPTV
               cubeHelper(x,y,z) = 1;
            end
         end
      end
      
   case {'spheric'}
      
      radiusPTV = xDim/16;
      
      for x = 1:xDim
         for y = 1:yDim
            for z = 1:zDim
               currPost = [x y z] - round([ct.cubeDim./2]);
               if  sqrt(sum(currPost.^2)) < radiusPTV
                  cubeHelper(x,y,z) = 1;
               end
            end
         end
      end
      
end



% extract the voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper);

% now we have ct data and cst data for a new phantom
display(ct);
display(cst);

%% Assign relative electron densities
vIxOAR = cst{ixOAR,4}{1};
vIxPTV = cst{ixPTV,4}{1};

ct.cubeHU{1}(vIxOAR) = 1;  % assign HU of water
ct.cubeHU{1}(vIxPTV) = 1;  % assign HU of water

%load(['patients' filesep 'BOXPHANTOM_TINY.mat'])
%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.
%%
% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon. In this
% example we would like to use photons for treatment planning. Next, we
% need to define a treatment machine to correctly load the corresponding 
% base data. matRad features generic base data in the file
% 'photons_Generic.mat'; consequently the machine has to be set to 'Generic'
pln.radiationMode = 'carbon';            
pln.machine       = 'GenericAPM';

%%
% Define the biological optimization model for treatment planning along
% with the quantity that should be used for optimization. Possible model values 
% are:
% 'none':     physical optimization;
% 'constRBE': constant RBE of 1.1; 
% 'MCN':      McNamara-variable RBE model for protons; 
% 'WED':      Wedenberg-variable RBE model for protons
% 'LEM':      Local Effect Model 
% and possible quantityOpt are 'physicalDose', 'effect' or 'RBExD'.
modelName    = 'LEM';
quantityOpt  = 'effect';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 20;
pln.propStf.gantryAngles  = [0 90];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% global INPUT_UNCERTAINTY;
% INPUT_UNCERTAINTY = 'phys'; %'phys', 'bio', 'biophys'
   
% retrieve bio model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

scenGenType  = 'apmScen';    % scenario creation type 'nomScen'  'wcScen' 'impScen' 'rndScen' 'apmScen'
% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType); 

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation
[cst,dij] = matRad_calcDose(ct,stf,pln,cst);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

 %% calculate variance
[cst,resultGUI] = matRad_calcVar(ct,cst,stf,pln,dij,resultGUI);
   
%%
cst{1,6}.robustness = 'PROB';
cst{2,6}.robustness = 'PROB';
for i= 1:size(cst,1)
    cst{i,4}{1} = cstOrg{i,4}{1};
end
param.w       = resultGUI.w;
resultGUIrob  = matRad_fluenceOptimization(dij,cst,pln,param);

%% calculate variance of robust pencil beam weights
[~,resultGUIrob] = matRad_calcVar(ct,cst,stf,pln,dij,resultGUIrob);

%% plot everthing
plane         = 3;
slice         = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
doseWindowExp = [0 max([max(max(resultGUI.RBExDExp(:,:,slice))) max(max(resultGUIrob.RBExDExpRob(:,:,slice)))])]*1.05;
doseWindowStd = [0 max([max(max(resultGUI.RBExDStdSingleFrac(:,:,slice))) max(max(resultGUIrob.RBExDStdSingleFracRob(:,:,slice)))])]*1.05;

figure,title('phantom plan')
subplot(221),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDExp,plane,slice,[],[],colorcube,[],doseWindowExp,[]);
subplot(222),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDStdSingleFrac,plane,slice,[],[],colorcube,[],doseWindowStd,[]);
subplot(223),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrob.RBExDExpRob,plane,slice,[],[],colorcube,[],doseWindowExp,[]);
subplot(224),matRad_plotSliceWrapper(gca,ct,cst,1,resultGUIrob.RBExDStdSingleFracRob,plane,slice,[],[],colorcube,[],doseWindowStd,[]);

%% append results and show them in GUI
resultGUI  = matRad_appendResultGUI(resultGUI,resultGUIrob,0);
matRadGUI

