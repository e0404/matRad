%% Test neutron dose calculation using PBKs - modification of example 1
clear
%% set matRad runtime configuration
matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%% Create a CT image series
xDim = 200;
yDim = 200;
zDim = 100;

ct.cubeDim      = [yDim xDim zDim]; % second cube dimension represents the x-coordinate
ct.resolution.x = 2;
ct.resolution.y = 2;
ct.resolution.z = 3;
ct.numOfCtScen  = 1;
 
% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;

ct.dicomInfo.RescaleIntercept = 1000;
ct.dicomInfo.RescaleSlope = 1;
%% Create the VOI data for the phantom
% Now we define structures a contour for the phantom and a target

ixOAR = 1;
ixPTV = 2;
ixBone = 3;
ixLung = 4;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'Body';
cst{ixOAR,3} = 'OAR';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'PTV_BNCT';
cst{ixPTV,3} = 'TARGET';
 
cst{ixBone,1} = 2;
cst{ixBone,2} = 'Bone';
cst{ixBone,3} = 'OAR';

cst{ixLung,1} = 3;
cst{ixLung,2} = 'Lung';
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
disp(ct);
disp(cst);


%% Assign relative electron densities
vIxOAR = cst{ixOAR,4}{1};
vIxPTV = cst{ixPTV,4}{1};
vIxBone = cst{ixBone,4}{1};
vIxLung = cst{ixLung,4}{1};

ct.cubeHU{1}(vIxOAR) = 0; % assign HU of water
ct.cubeHU{1}(vIxPTV) = 0; % assign HU of water
ct.cubeHU{1}(vIxBone) = 2000; % assign HU of water
ct.cubeHU{1}(vIxLung) = -300; % assign HU of water

% Clean up
clearvars -except ct cst

%% Treatment Plan
pln.radiationMode = 'neutrons';
pln.machine       = 'BNCT';
pln.propDoseCalc.engine = 'MCNP';       

pln.numOfFractions        = 1;
pln.propStf.gantryAngles  = 0; %[0 20 340];
pln.propStf.couchAngles   = 0; %[0 0 0];
pln.propStf.bixelWidth    = 50;

% Dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 2; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% Set add margin to false in order to avoid an oversized margin leading to
% large number of rays
pln.propStf.addMargin = false;

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);
% cst{2,6}{1,1}.parameters{1,1} = 1.5;

%% Dose Calculation
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for intensity-modulated photon therapy
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the
% clinical objectives and constraints underlying the radiation treatment.
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI