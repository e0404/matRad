% configuration script for uncertainty sampling
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% examine only specific structures? improves calculation time
examineStructures = {}; % e.g. examinedStructures = {'CTV', 'OAR'};

multScen = matRad_multScen([],'impScen');
% a) define shift scenarios
multScen.numOfShiftScen       = [0 0 0];          % number of shifts in x y and z direction       
multScen.shiftSize            = [4.5 4.5 4.5];          % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
multScen.shiftGenType         = 'equidistant';    % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftCombType        = 'individual';     % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen)
                                                  % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
                                                  % combined: number of shifts in each directions must be the same. only use when shifts are sampled

                                             
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

%% report parameters
param.operator = 'matrad';

% default set of percentiles for scenario analysis; uncomment to override
% param.percentiles = [0.01 0.05 0.125 0.25 0.5 0.75 0.875 0.95 0.99];
% default criteria for gamma analysis between nominal and mean dose
% param.criteria = [3 3]; %%% [X % Y mm]

%% start calculation
matRad_calcStudy(examineStructures, multScen, [], param);

% exit;
