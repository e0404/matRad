%% Example: Brachytheraphy Treatment Plan
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% List of contents
% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup an HDR brachy dose calculation and 
% (iii) how to inversely optimize holding position intensties
% (iv) how to visually and quantitatively evaluate the result
% (v) how to verify that functions do the right thing

%% I Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_rc;
load 'PROSTATE.mat';

%% I - update/set dose objectives for brachytherapy
% The sixth column represents dose objectives and constraints for
% optimization: First, the objective function for the individual structure
% is chosen, its parameters denote doses that should not be tranceeded
% towards higher or lower doses (SquaredOverdose, SquaredUnderdose) or
% doses that are particularly aimed for (SquaredUnderDose).

disp(cst{6,6}{1});

% Following frequently prescribed planning doses of 15 Gy
% (https://pubmed.ncbi.nlm.nih.gov/22559663/) objectives can be updated to:

% the planning target was changed to the clinical segmentation of the
% prostate bed.
cst{5,3}    = 'TARGET';
cst{5,6}{1} = struct(DoseObjectives.matRad_SquaredUnderdosing(100,15));
cst{5,6}{2} = struct(DoseObjectives.matRad_SquaredOverdosing(100,17.5));
cst{6,5}.Priority = 1;

% In this example, the lymph nodes will not be part of the treatment:
cst{7,6}    =  [];
cst{7,3}    =  'OAR';

%A PTV is not needed, but we will use it to control the dose fall off
cst{6,3}    =  'OAR';
cst{6,6}{1} =  struct(DoseObjectives.matRad_SquaredOverdosing(100,12));
cst{6,5}.Priority = 2;

% Rectum Objective
cst{1,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));

% Bladder Objective
cst{8,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,7.5));

% Body NT objective
cst{9,6}{1} = struct(DoseObjectives.matRad_MeanDose(1));


%% II.1 Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

% First of all, we need to define what kind of radiation modality we would
% like to use. Possible values are photons, protons or carbon 
% (external beam) or brachy as an invasive tratment option.
% In this case we want to use brachytherapy. Then, we need to define a 
% treatment machine to correctly load the corresponding base data.
% matRad includes example base data for HDR and LDR brachytherapy.
% Here we will use HDR. By this means matRad will look for 'brachy_HDR.mat'
% in our root directory and will use the data provided in there for 
% dose calculation.

pln.radiationMode   = 'brachy'; 
pln.machine         = 'HDR';    % 'LDR' or 'HDR' for brachy


%% II.1 - needle and template geometry
% Now we have to set some parameters for the template and the needles. 
% Let's start with the needles: Seed distance is the distance between
% two neighbouring seeds or holding points on one needle or catheter. The
% seeds No denotes how many seeds/holding points there are per needle.

pln.propStf.needle.seedDistance      = 10; % [mm] 
pln.propStf.needle.seedsNo           = 6; 

%% II.1 - template position
% The implantation is normally done through a 13 x 13 template from the 
% patients inferior, which is the negative z axis here.
% The direction of the needles is defined by template normal.
% Neighbour distances are called by bixelWidth, because this field is also
% used for external beam therapy.
% The needles will be positioned right under the target volume pointing up.

pln.propStf.template.normal      = [0,0,1];
pln.propStf.bixelWidth   = 5; % [mm] template grid distance
pln.propStf.templateRoot = matRad_getTemplateRoot(ct,cst); % mass center of
% target in x and y and bottom in z

% Here, we define active needles as 1 and inactive needles
% as 0. This is the x-y plane and needles point in z direction. 
% A checkerboard pattern is frequantly used. The whole geometry will become
% clearer when it is displayed in 3D view in the next section.

pln.propStf.template.activeNeedles = [0 0 0 1 0 1 0 1 0 1 0 0 0;... % 7.0
                                      0 0 1 0 1 0 0 0 1 0 1 0 0;... % 6.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 6.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 5.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 5.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 4.0
                                      1 0 1 0 1 0 0 0 1 0 1 0 1;... % 4.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 3.0
                                      1 0 1 0 1 0 1 0 1 0 1 0 1;... % 2.5
                                      0 1 0 1 0 1 0 1 0 1 0 1 0;... % 2.0
                                      1 0 1 0 1 0 0 0 0 0 1 0 1;... % 1.5
                                      0 0 0 0 0 0 0 0 0 0 0 0 0];   % 1.0
                                     %A a B b C c D d E e F f G

pln.propStf.isoCenter    = matRad_getIsoCenter(cst,ct,0); %  target center



%% II.1 - dose calculation options
% for dose calculation we use eather the 2D or the 1D formalism proposed by
% TG 43. Also, set resolution of dose calculation and optimization.
% If your system gets stuck with the resolution, you can lower it to 10 or
% 20, just to get an initial result. Otherwise, reduce the number of
% needles.
% Calculation time will be reduced by one tenth when we define a dose
% cutoff distance.


pln.propDoseCalc.TG43approximation = '2D'; %'1D' or '2D' 

pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]

pln.propDoseCalc.DistanceCutoff    = 130; %[mm] sets the maximum distance
                                            %to which dose is calculated. 

% the standard interior point optimizer IPOPT can be used
pln.propOpt.optimizer       = 'IPOPT';

%% II.1 - book keeping
% Some field names have to be kept although they don't have a direct
% relevance for brachy therapy.
pln.propOpt.bioOptimization = 'none';
pln.propOpt.runDAO          = false;  
pln.propOpt.runSequencing   = false; 
pln.propStf.gantryAngles    = []; 
pln.propStf.couchAngles     = []; 
pln.propStf.numOfBeams      = 0;
pln.numOfFractions          = 1; 

%% II.1 - view plan
% Et voila! Our treatment plan structure is ready. Lets have a look:
disp(pln);


%% II.2 Steering Seed Positions From STF
% The steering file struct contains all needls/catheter geometry with the
% target volume, number of needles, seeds and the positions of all needles
% The one in the end enables visualization.

stf = matRad_generateStf(ct,cst,pln,1);

%% II.2 - view stf
% The 3D view is interesting, but we also want to know how the stf struct
% looks like.

disp(stf)

%% II.3 - Dose Calculation
% Let's generate dosimetric information by pre-computing a dose influence 
% matrix for seed/holding point intensities. Having dose influences
% available allows subsequent inverse optimization.
% Don't get inpatient, this can take a few seconds...

dij = matRad_calcBrachyDose(ct,stf,pln,cst);

%% III Inverse Optimization for brachy therapy
% The goal of the fluence optimization is to find a set of holding point
% times which yield the best possible dose distribution according to
% the clinical objectives and constraints underlying the radiation 
% treatment. Once the optimization has finished, trigger to 
% visualize the optimized dose cubes.

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
matRadGUI;

%% IV.1 Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice

slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet);

%% IV.2 Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
[dvh,qi]               = matRad_indicatorWrapper(cst,pln,resultGUI);

%% V Unit testing
% It is very useful to test individual function for against results that
% have been calculated using a trusted other method.
% run the unit test suites for geometry and dose calculation

matRad_runBrachyTestSuite

% if you want to built the code further, it is recommended, to write at
% least one unit test for each complicated function.