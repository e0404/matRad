%% Test neutron dose calculation using PBKs - modification of example 1
clear
%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Create a CT image series
xDim = 200;
yDim = 200;
zDim = 50;

ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate
ct.resolution.x = 2;
ct.resolution.y = 2;
ct.resolution.z = 3;
ct.numOfCtScen  = 1;
 
% create a ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;

%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target

ixOAR = 1;
ixPTV = 2;
ixBone = 3;
ixLung = 4;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';
 
cst{ixBone,1} = 2;
cst{ixBone,2} = 'bone';
cst{ixBone,3} = 'OAR';

cst{ixLung,1} = 3;
cst{ixLung,2} = 'lung';
cst{ixLung,3} = 'OAR';

% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass  = 1;
cst{ixOAR,5}.alphaX       = 0.1000;
cst{ixOAR,5}.betaX        = 0.0500;
cst{ixOAR,5}.Priority     = 2;
cst{ixOAR,5}.Visible      = 1;
cst{ixOAR,5}.visibleColor = [0 0 0];

cst{ixBone,5}.TissueClass  = 1;
cst{ixBone,5}.alphaX       = 0.1000;
cst{ixBone,5}.betaX        = 0.0500;
cst{ixBone,5}.Priority     = 2;
cst{ixBone,5}.Visible      = 1;
cst{ixBone,5}.visibleColor = [1 0 0];

cst{ixLung,5}.TissueClass  = 1;
cst{ixLung,5}.alphaX       = 0.1000;
cst{ixLung,5}.betaX        = 0.0500;
cst{ixLung,5}.Priority     = 2;
cst{ixLung,5}.Visible      = 1;
cst{ixLung,5}.visibleColor = [0 1 0];

% define objective as struct for compatibility with GNU Octave I/O
cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));
cst{ixBone,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));
cst{ixLung,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));

cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,5}.visibleColor = [1 1 1];

% define objective as struct for compatibility with GNU Octave I/O
cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));

%% Lets create either a cubic or a spheric phantom
cubeHelper = zeros(ct.cubeDim);

% Soft tissue OAR
      xLowOAR  = round(xDim/2 - xDim/4);
      xHighOAR = round(xDim/2 + xDim/4);
      yLowOAR  = round(yDim/2 - yDim/4);
      yHighOAR = round(yDim/2 + yDim/4);
      zLowOAR  = round(zDim/2 - zDim/4);
      zHighOAR = round(zDim/2 + zDim/4);
      
      for x = xLowOAR:1:xHighOAR
         for y = yLowOAR:1:yHighOAR
            for z = zLowOAR:1:zHighOAR
               cubeHelper(y,x,z) = 1;
            end
         end
      end

% radiusOAR = xDim/4;
% 
% for x = 1:xDim
%     for y = 1:yDim
%         for z = 1:zDim
%             currPost = [y x z] - round([ct.cubeDim./2]);
%             if  sqrt(sum(currPost.^2)) < radiusOAR
%                 cubeHelper(y,x,z) = 1;
%             end
%         end
%     end
% end
clear currPost

% extract the voxel indices and save it in the cst
cst{ixOAR,4}{1} = find(cubeHelper);

% Bone
cubeHelper = zeros(ct.cubeDim);

radiusBone = xDim/12;

for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [y x z] - round([ct.cubeDim./2]);
            currPost(2) = currPost(2) + round([ct.cubeDim(2)./8]);
            if  sqrt(sum(currPost.^2)) < radiusBone
                cubeHelper(y,x,z) = 1;
            end
        end
    end
end
clear currPost

% extract the voxel indices and save it in the cst
cst{ixBone,4}{1} = find(cubeHelper);

% Lung
cubeHelper = zeros(ct.cubeDim);

radiusLung = xDim/12;

for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [y x z] - round([ct.cubeDim./2]);
            currPost(2) = currPost(2) - round([ct.cubeDim(2)./8]);
            if  sqrt(sum(currPost.^2)) < radiusLung
                cubeHelper(y,x,z) = 1;
            end
        end
    end
end
clear currPost
% extract the voxel indices and save it in the cst
cst{ixLung,4}{1} = find(cubeHelper);

% PTV
cubeHelper = zeros(ct.cubeDim);
radiusPTV = xDim/24;

for x = 1:xDim
    for y = 1:yDim
        for z = 1:zDim
            currPost = [x y z] - round([ct.cubeDim./2]);
            currPost(2) = currPost(2) + round([ct.cubeDim(2)./5]);
            if  sqrt(sum(currPost.^2)) < radiusPTV
                cubeHelper(y,x,z) = 1;
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
vIxBone = cst{ixBone,4}{1};
vIxLung = cst{ixLung,4}{1};

ct.cubeHU{1}(vIxOAR) = 0; % assign HU of water
ct.cubeHU{1}(vIxPTV) = 0; % assign HU of water
ct.cubeHU{1}(vIxBone) = 2000; % assign HU of water
ct.cubeHU{1}(vIxLung) = -300; % assign HU of water

ct.dicomInfo.ManufacturerModelName = 'neutronXScorrfinalCalcMLC1';
ct.dicomInfo.ConvolutionKernel = 'neutronField';
ct.dicomInfo.Manufacturer = 'matRad';

ct.dicomInfo_gammaComp.ManufacturerModelName = 'allWater';
ct.dicomInfo_gammaComp.ConvolutionKernel = 'neutronField';
ct.dicomInfo_gammaComp.Manufacturer = 'matRad';

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
pln.radiationMode = 'neutrons'; %'photons'; %'neutrons';            
pln.machine       = 'mixedField_MLCv1_MEDAPP'; %'Generic'; %'neutronField_FRM2_MLCv1';

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
modelName    = 'constRBE'; %'constRBE'; %'none';
quantityOpt  = 'RBExD'; %'RBExD'; %'physicalDose';                                             

%%
% The remaining plan parameters are set like in the previous example files
pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = [20]; %; 20];
pln.propStf.couchAngles   = [0]; %; 10];
pln.propStf.bixelWidth    = 20; %90; %5; %90;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% retrieve bio model parameters
pln.bioParam.model = 'constRBE';
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt,modelName);

% retrieve nominal scenario for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,'nomScen'); 

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

pln.propDoseCalc.omitFNTmixedField = false;
pln.propDoseCalc.useCustomPrimaryNeutronFluence = false; %true;
pln.propDoseCalc.geometricCutOff = 300;

pln.propDoseCalc.doseCorrectionKERMAfactors = 'correctionFactorsKERMA-MEDAPP';
pln.propDoseCalc.neutronOverGammaRatioRefDepth = 1.82; % neutron over gamma ratio in 5 cm depth MLC 3

% Set add margin to false in order to avoid an oversized margin leading to
% large number of rays
pln.propStf.addMargin = false;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
cst{2,6}{1,1}.parameters{1,1} = 1.5;

%% Dose Calculation
dij = matRad_calcNeutronDose(ct,stf,pln,cst);

%% Export dij matrix
%matRad_exportDij('dij.bin',dij,stf);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI

% %% Plot the resulting dose slice
% plane      = 3;
% slice      = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
% doseWindow = [0 max([resultGUI.physicalDose(:)])];
% 
% figure,title('phantom plan')
% matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],[],colorcube,[],doseWindow,[]);
% 
% 
% %% 
% % We export the the created phantom & dose as dicom. This is handled by the 
% % class matRad_DicomExporter. When no arguments are given, the exporter searches
% % the workspace itself for matRad-structures. The output directory can be set by
% % the property dicomDir. While the different DICOM datasets (ct, RTStruct, etc) 
% % can be exported individually, we call the wrapper to do all possible exports.
% dcmExport = matRad_DicomExporter();
% dcmExport.dicomDir = [pwd filesep 'dicomExport'];
% dcmExport.matRad_exportDicom();

