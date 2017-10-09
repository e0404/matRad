%% script to generate trusted scenarios
warning('You need to recheck everything by hand before generating a trusted testCase file. Copy hash to testScen.');

ct.numOfCtScen = 1;
pln.sampling = true;
i = 1;
%%
% a) define shift scenarios
multScen.numOfShiftScen       = [0 0 0];          % number of shifts in x y and z direction       
multScen.shiftSize            = [3 3 3];          % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
multScen.shiftGenType         = 'equidistant';    % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
multScen.shiftCombType        = 'individual';     % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
                                                  % combined:    combine shift scenarios;                 number of shift scenarios is multScen.numOfShiftScen(1)
                                                  % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
multScen.shiftGen1DIsotropy   = '+-';             % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

% b) define range error scenarios                                                
multScen.numOfRangeShiftScen  = 4; % number of absolute and/or relative range scnearios. 
                                                  % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
multScen.maxAbsRangeShift     = 0;                % maximum absolute over and undershoot in mm   
multScen.maxRelRangeShift     = 7;              % maximum relative over and undershoot in % 
multScen.rangeCombType        = 'combined';       % individual: no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
multScen.rangeGenType         = 'equidistant';    % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
multScen.scenCombType         = 'individual';     % individual:  no combination of range and setup scenarios, 
                                                  % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                  % permuted:    create every possible combination of range and setup scenarios
multScen.includeNomScen       = false;

% define standard deviation of normal distribution - important for probabilistic treatment planning
multScen.rangeRelSD           = 3.5;               % given in [%]   
multScen.rangeAbsSD           = 1;                 % given in [mm]   
multScen.shiftSD              = [2 2 2];           % given in [mm]

multScenCases(i).multScen = multScen;
pln = matRad_setPlanUncertainties(ct, pln, multScen);
multScenCases(i).prob = pln.multScen.scenProb;
multScenCases(i).scenParam = pln.multScen.scenForProb;

clear multScen;
i = i + 1;

%% write to file
testCase.multScenCases = multScenCases;
DataHash(testCase)
save(fullfile('unitTest','testCases.mat'), 'testCase');