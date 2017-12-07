%% sampling
% examine only specific structures? improves calculation time
examineStructures = {}; % e.g. examinedStructures = {'CTV', 'OAR'};

% a) define shift scenarios
multScen.shiftSize            = 2.5;                % maximum shift in multiples of shiftSD (grid)
multScen.shiftGenType         = 'grid';         % grid: equidistant shifts, sampled: sample shifts from multivariate normal distribution
multScen.numOfShiftScen       = 729;              % number of shifts for grid use (64,125,216,343,512,729,1000,N^3)


% b) define range error scenarios                                                
multScen.numOfRangeShiftScen  = 0; % number of absolute and/or relative range scnearios. 
                                                  % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
multScen.maxAbsRangeShift     = 0;                % maximum absolute over and undershoot in mm   
multScen.maxRelRangeShift     = 7;                % maximum relative over and undershoot in % 
multScen.rangeCombType        = 'combined';       % serial: no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
multScen.rangeGenType         = 'equidistant';    % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution

multScen.scenCombType         = 'individual';     % serial:  no combination of range and setup scenarios, 
                                                  % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                  % permuted:    create every possible combination of range and setup scenarios
multScen.includeNomScen       = false;

% define standard deviation of normal distribution - important for probabilistic treatment planning
multScen.shiftSD              = [1 100 10000];     % given in [mm]
multScen.rangeRelSD           = 3.5;               % given in [%]   
multScen.rangeAbsSD           = 1;                 % given in [mm]   

%% path for output pdf and mat
param.outputPath = pwd;
% addpath(genpath('C:\git\matRad')) % optional add your matRad path here if not yet added to searchpath

%% add your name here
param.operator = 'Werner Heisenberg';

%% start calculation
matRad_calcStudy(examineStructures, multScen, param);

% exit;
