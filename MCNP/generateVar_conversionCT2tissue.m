% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decription: Generate variable 'conversionCT2tissue.mat' where HU
% intervals for tissue identification are defined
%
% References: 
%   [1] Schneider  et al. Correlation between CT numbers and
%       tissue parameters needed for Monte Carlo simulations of clinical 
%       dose distribution, Phys. Med. Biol. 45 (2000)
%   [2] Schneider et al. The calibration of CT HU for RT treatment
%       planning, Phys. Med. Bio. 41 (1996)
%   [3] ICRU 46
%   [4] MCNP6 Manual Part 3 Appendix G
%   [5] DeMarco et al. A CTbased Monte Carlo simulation tool for dosimetry
%       planning and analysis. Medical physics 25.1 (1998)
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 10/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;


binIntervals(1).name = 'air';
binIntervals(2).name = 'lung';
binIntervals(3).name = 'softTissue';
binIntervals(4).name = 'bone';
binIntervals(5).name = 'skin';
binIntervals(6).name = 'ptv_bnct';

% Define parameter for segementation
segVar.upperLimitAir = 50; % upper limit for scaled HU below which all material is assigned to air, chosen according to Schneider et al. (2000)
segVar.densityAir = 1.205e-3; % density air in g/cm^3 from NIST database

binIntervals(1).HUbin = [0 segVar.upperLimitAir];
binIntervals(2).HUbin = [segVar.upperLimitAir 830];     % lung tissue has HU_lung = 259 -> see paper Schneider et al. (1996)
binIntervals(3).HUbin = [830 1280];                     % all limits according to Schneider et al. 1996/2000 except upper limit soft tissue according to DeMarco et al. 1998
binIntervals(4).HUbin = [1280 5000];
binIntervals(5).HUbin = [];
binIntervals(6).HUbin = [];

binIntervals(1).importance = 1; % importance for particles: 'MODE N P E H D T S A #' is automatically adjusted for KERMA calculations
binIntervals(2).importance = 5; 
binIntervals(3).importance = 5;
binIntervals(4).importance = 5;
binIntervals(5).importance = 5;
binIntervals(6).importance = 5;


%% Define the compositon of the material used for the simulation according to ICRU Report 46 and MCNP Appendix G for cross sections

%% Cross sections

crossSectionsLibrary(1,:)='.66c';
crossSectionsLibrary(2,:)= '.80c'; 

%% Material Definition: Isotopes only added, when rounded percentage is at least 1%
%% Air - Taken from Schneider et al. (2000)
binIntervals(1).ZAID(1,1) = 7014; % N 
binIntervals(1).crossSection(1,1) = 2;
binIntervals(1).percentageMass(1,1) = -0.755;   % 

binIntervals(1).ZAID(1,2) = 8016; % O
binIntervals(1).crossSection(1,2) = 2;
binIntervals(1).percentageMass(1,2) = -0.232;

binIntervals(1).ZAID(1,3) = 18040; % Ar
binIntervals(1).crossSection(1,3) = 2;
binIntervals(1).percentageMass(1,3) = -0.013;

%% Lung 
binIntervals(2).ZAID(1,1) = 1001; % H from ICRU Report 46
binIntervals(2).crossSection(1,1) = 1;
binIntervals(2).percentageMass(1,1) = -0.103;

binIntervals(2).ZAID(1,2) = 6000; % C from ICRU Report 46
binIntervals(2).crossSection(1,2) = 2;
binIntervals(2).percentageMass(1,2) = -0.105;

binIntervals(2).ZAID(1,3) = 7014; % N from ICRU Report 46
binIntervals(2).crossSection(1,3) = 2;
binIntervals(2).percentageMass(1,3) = -0.031;

binIntervals(2).ZAID(1,4) = 8016; % O from ICRU Report 46
binIntervals(2).crossSection(1,4) = 2;
binIntervals(2).percentageMass(1,4) = -0.749;

binIntervals(2).ZAID(1,5) = 11023; % Na from ICRU Report 46
binIntervals(2).crossSection(1,5) = 2;
binIntervals(2).percentageMass(1,5) = -0.002;

binIntervals(2).ZAID(1,6) = 15031; % P from ICRU Report 46
binIntervals(2).crossSection(1,6) = 2;
binIntervals(2).percentageMass(1,6) = -0.002;

binIntervals(2).ZAID(1,7) = 16032; % S from ICRU Report 46
binIntervals(2).crossSection(1,7) = 2;
binIntervals(2).percentageMass(1,7) = -0.003*0.95;
binIntervals(2).ZAID(1,8) = 16033; % S from ICRU Report 46
binIntervals(2).crossSection(1,8) = 2;
binIntervals(2).percentageMass(1,8) = -0.003*0.01;
binIntervals(2).ZAID(1,9) = 16034; % S from ICRU Report 46
binIntervals(2).crossSection(1,9) = 2;
binIntervals(2).percentageMass(1,9) = -0.003*0.04;

binIntervals(2).ZAID(1,10) = 17035; % Cl from ICRU Report 46
binIntervals(2).crossSection(1,10) = 2;
binIntervals(2).percentageMass(1,10) = -0.003*0.76;
binIntervals(2).ZAID(1,11) = 17037; % Cl from ICRU Report 46
binIntervals(2).crossSection(1,11) = 2;
binIntervals(2).percentageMass(1,11) = -0.003*0.24;

binIntervals(2).ZAID(1,12) = 19039; % K from ICRU Report 46
binIntervals(2).crossSection(1,12) = 2;
binIntervals(2).percentageMass(1,12) = -0.002*0.93;
binIntervals(2).ZAID(1,13) = 19041; % K from ICRU Report 46
binIntervals(2).crossSection(1,13) = 2;
binIntervals(2).percentageMass(1,13) = -0.002*0.07;

%% Soft Tissue
binIntervals(3).ZAID(1,1) = 1001; % H from ICRU Report 46
binIntervals(3).crossSection(1,1) = 1;
binIntervals(3).percentageMass(1,1) = -0.101;

binIntervals(3).ZAID(1,2) = 6000; % C from ICRU Report 46
binIntervals(3).crossSection(1,2) = 2;
binIntervals(3).percentageMass(1,2) = -0.111;

binIntervals(3).ZAID(1,3) = 7014; % N from ICRU Report 46
binIntervals(3).crossSection(1,3) = 2;
binIntervals(3).percentageMass(1,3) = -0.026;

binIntervals(3).ZAID(1,4) = 8016; % O from ICRU Report 46
binIntervals(3).crossSection(1,4) = 2;
binIntervals(3).percentageMass(1,4) = -0.762;

%% Bone
binIntervals(4).ZAID(1,1) = 1001; % H from ICRU Report 44; all taken from DeMarco et al. 1998
binIntervals(4).crossSection(1,1) = 1;
binIntervals(4).percentageMass(1,1) = -0.034;

binIntervals(4).ZAID(1,2) = 6000; % C from ICRU Report 46
binIntervals(4).crossSection(1,2) = 2;
binIntervals(4).percentageMass(1,2) = -0.155;

binIntervals(4).ZAID(1,3) = 7014; % N from ICRU Report 46
binIntervals(4).crossSection(1,3) = 2;
binIntervals(4).percentageMass(1,3) = -0.042;

binIntervals(4).ZAID(1,4) = 8016; % O from ICRU Report 44
binIntervals(4).crossSection(1,4) = 2;
binIntervals(4).percentageMass(1,4) = -0.435;

binIntervals(4).ZAID(1,5) = 15031; % P from ICRU Report 44
binIntervals(4).crossSection(1,5) = 2;
binIntervals(4).percentageMass(1,5) = -0.103;

binIntervals(4).ZAID(1,6) = 20040; % Ca from ICRU Report 44
binIntervals(4).crossSection(1,6) = 2;
binIntervals(4).percentageMass(1,6) = -0.225*0.97;
binIntervals(4).ZAID(1,7) = 20042; % Ca from ICRU Report 44
binIntervals(4).crossSection(1,7) = 2;
binIntervals(4).percentageMass(1,7) = -0.225*0.01;
binIntervals(4).ZAID(1,8) = 20044; % Ca from ICRU Report 44
binIntervals(4).crossSection(1,8) = 2;
binIntervals(4).percentageMass(1,8) = -0.225*0.02;

%% Skin 
binIntervals(5).ZAID(1,1) = 1001; % H from ICRU Report 46
binIntervals(5).crossSection(1,1) = 1;
binIntervals(5).percentageMass(1,1) = -0.100;

binIntervals(5).ZAID(1,2) = 6000; % C from ICRU Report 46
binIntervals(5).crossSection(1,2) = 2;
binIntervals(5).percentageMass(1,2) = -0.204*0.99;

binIntervals(5).ZAID(1,3) = 7014; % N from ICRU Report 46
binIntervals(5).crossSection(1,3) = 2;
binIntervals(5).percentageMass(1,3) = -0.042;

binIntervals(5).ZAID(1,4) = 8016; % O from ICRU Report 46
binIntervals(5).crossSection(1,4) = 2;
binIntervals(5).percentageMass(1,4) = -0.645;

binIntervals(5).ZAID(1,5) = 11023; % Na from ICRU Report 46
binIntervals(5).crossSection(1,5) = 2;
binIntervals(5).percentageMass(1,5) = -0.002;

binIntervals(5).ZAID(1,6) = 15031; % P from ICRU Report 44
binIntervals(5).crossSection(1,6) = 2;
binIntervals(5).percentageMass(1,6) = -0.001;

binIntervals(5).ZAID(1,7) = 16032; % S from ICRU Report 46
binIntervals(5).crossSection(1,7) = 2;
binIntervals(5).percentageMass(1,7) = -0.002*0.95;
binIntervals(5).ZAID(1,8) = 16033; % S from ICRU Report 46
binIntervals(5).crossSection(1,8) = 2;
binIntervals(5).percentageMass(1,8) = -0.002*0.01;
binIntervals(5).ZAID(1,9) = 16034; % S from ICRU Report 46
binIntervals(5).crossSection(1,9) = 2;
binIntervals(5).percentageMass(1,9) = -0.002*0.04;

binIntervals(5).ZAID(1,10) = 17035; % Cl from ICRU Report 46
binIntervals(5).crossSection(1,10) = 2;
binIntervals(5).percentageMass(1,10) = -0.003*0.76;
binIntervals(5).ZAID(1,11) = 17037; % Cl from ICRU Report 46
binIntervals(5).crossSection(1,11) = 2;
binIntervals(5).percentageMass(1,11) = -0.003*0.24;

binIntervals(5).ZAID(1,12) = 19039; % K from ICRU Report 46
binIntervals(5).crossSection(1,12) = 2;
binIntervals(5).percentageMass(1,12) = -0.001*0.93;
binIntervals(5).ZAID(1,13) = 19041; % K from ICRU Report 46
binIntervals(5).crossSection(1,13) = 2;
binIntervals(5).percentageMass(1,13) = -0.001*0.07;

%% Soft Tissue with Boron Content
% Fraction of mass necessary for BNCT (Chandra et al., 2015): 20mug/g

boronPercentPerMass = 30e-6; %30e-6; %20e-6;

binIntervals(6).ZAID(1,1) = 1001; % H from ICRU Report 46
binIntervals(6).crossSection(1,1) = 1;
binIntervals(6).percentageMass(1,1) = -0.101 * (1 - boronPercentPerMass);

binIntervals(6).ZAID(1,2) = 6000; % C from ICRU Report 46
binIntervals(6).crossSection(1,2) = 2;
binIntervals(6).percentageMass(1,2) = -0.105*0.99 * (1 - boronPercentPerMass);

binIntervals(6).ZAID(1,3) = 7014; % N from ICRU Report 46
binIntervals(6).crossSection(1,3) = 2;
binIntervals(6).percentageMass(1,3) = -0.026 * (1 - boronPercentPerMass);

binIntervals(6).ZAID(1,4) = 8016; % O from ICRU Report 46
binIntervals(6).crossSection(1,4) = 2;
binIntervals(6).percentageMass(1,4) = -0.762 * (1 - boronPercentPerMass);

binIntervals(6).ZAID(1,5) = 5010; % B-10
binIntervals(6).crossSection(1,5) = 2;
binIntervals(6).percentageMass(1,5) = -boronPercentPerMass;

%% Controle whether percentages are all one
dummySummy = 0;
for i=1:size(binIntervals,2)
    dummySummy = dummySummy + sum(binIntervals(i).percentageMass(1,:));
end
clear i;
if abs(dummySummy)~=size(binIntervals,2)
    warning('One of the material definitions is missing one or more components. Change here otherwise MCNP will normalize automatically...')
end
clear dummySummy;

%% Save
save('conversionCT2tissue.mat', 'binIntervals', 'segVar','crossSectionsLibrary')
