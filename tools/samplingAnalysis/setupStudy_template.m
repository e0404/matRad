%% sampling
% examine only specific structures? improves calculation time
examineStructures = {}; % e.g. examinedStructures = {'CTV', 'OAR'};

multScen = matRad_multScen([],'impScen');
% a) define shift scenarios
multScen.numOfShiftScen       = [0 0 0];          % number of shifts in x y and z direction       
multScen.shiftSize            = [4.5 4.5 4.5];          % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
multScen.shiftGenType         = 'equidistant';    % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftCombType        = 'individual';     % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen)
                                                  % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
% b) define range error scenarios                                                
multScen.numOfRangeShiftScen  = 30;                % number of absolute and/or relative range scnearios. 
                                                  % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
multScen.maxAbsRangeShift     = 2;                % maximum absolute over and undershoot in mm   
multScen.maxRelRangeShift     = 4.5;                % maximum relative over and undershoot in % 
multScen.rangeCombType        = 'combined';       % individual: no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
multScen.rangeGenType         = 'equidistant';    % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
multScen.scenCombType         = 'individual';     % individual:  no combination of range and setup scenarios, 
                                                  % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                  % permuted:    create every possible combination of range and setup scenarios
multScen.includeNomScen       = false;


% define standard deviation of normal distribution - important for probabilistic treatment planning
multScen.rangeRelSD           = 0;               % given in [%]   
multScen.rangeAbsSD           = 1.8;                 % given in [mm]   
multScen.shiftSD              = [1.8 1.8 1.8];           % given in [mm]

multScen = multScen.matRad_createValidInstance();

%% path for output pdf and mat
param.outputPath = pwd;
% addpath(genpath('C:\git\matRad')) % optional add your matRad path here if not yet added to searchpath

%% report parameters
param.operator = 'Lucas-Raphael Mueller';

% default set of percentiles for scenario analysis; uncomment to override
% param.percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
% default criteria for gamma analysis between nominal and mean dose
% param.criteria = [3 3]; %%% [X % Y mm]

%% start calculation
matRad_calcStudy(examineStructures, multScen, [], param);

% exit;
