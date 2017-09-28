function [pln] = matRad_setPlanUncertainties(ct,pln, multScen, param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_setPlanUncertainties function provides functionalities to define 
% treatment planning uncertainties
% 
% call
%   [cst,pln] = matRad_setPlanUncertainties(ct,cst,pln)
%
% input
%   ct:             ct cube
%   pln:            matRad plan meta information struct
%
% output
%   pln:            matRad's plan meta information struct including a sub-structure 
%                   pln.multScen holding information about multiple treatment plan scenarios
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


if ~isfield(pln,'robOpt')
   pln.robOpt = false;
end

if ~isfield(pln,'sampling')
   pln.sampling = false;
end

% define standard deviation of normal distribution - only relevant for probabilistic treatment planning
requiredFields = {'rangeRelSD', 'rangeAbsSD', 'shiftSD'};
if ~exist('multScen', 'var') || sum(isfield(multScen,requiredFields)) == 0
    multScen.rangeRelSD           = 3.5;               % given in %
    multScen.rangeAbsSD           = 1;                 % given in [mm]   
    multScen.shiftSD              = [3 3 3];           % given in [mm]
end


 %% create multiple scenario struc

% define parameters for individual treatment planning scenarios 
multScen.numOfCtScen          = ct.numOfCtScen; % number of imported ct scenarios


% if robust optimization is disabled then create solely the nominal scenario
if ~pln.robOpt && ~pln.sampling
    % a) define shift scenarios
   multScen.numOfShiftScen       = [0 0 0];         % number of shifts in x y and z direction       
   multScen.shiftSize            = [0 0 0];         % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
   multScen.shiftGenType         = 'equidistant';   % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.shiftCombType        = 'individual';    % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
                                                    % combined:    combine shift scenarios;                 number of shift scenarios is multScen.numOfShiftScen(1)
                                                    % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
   multScen.shiftGen1DIsotropy   = '+-';            % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

   % b) define range error scenarios                                                
   multScen.numOfRangeShiftScen  = 0;               % number of absolute and/or relative range scnearios. 
                                                    % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
   multScen.maxAbsRangeShift     = 0;               % maximum absolute over and undershoot in mm   
   multScen.maxRelRangeShift     = 0  ;             % maximum relative over and undershoot in % 
   multScen.rangeCombType        = 'combined';      % individual: no combination of absolute and relative range scenarios
                                                    % combined:    combine absolute and relative range scenarios
   multScen.rangeGenType         = 'equidistant';   % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
   multScen.scenCombType         = 'individual';    % individual:  no combination of range and setup scenarios, 
                                                    % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                    % permuted:    create every possible combination of range and setup scenarios
   multScen.includeNomScen       = true;
   
% definition of scenarios for robust optimization   
elseif pln.robOpt && ~pln.sampling 
   
   % a) define shift scenarios
   multScen.numOfShiftScen       = [2 2 2];         % number of shifts in x y and z direction       
   multScen.shiftSize            = [3 3 3];         % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
   multScen.shiftGenType         = 'equidistant';   % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.shiftCombType        = 'individual';    % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
                                                    % combined:    combine shift scenarios;                 number of shift scenarios is multScen.numOfShiftScen(1)
                                                    % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
   multScen.shiftGen1DIsotropy   = '+-';            % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

   % b) define range error scenarios                                                
   multScen.numOfRangeShiftScen  = 2;               % number of absolute and/or relative range scnearios. 
                                                    % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
   multScen.maxAbsRangeShift     = 1;               % maximum absolute over and undershoot in mm   
   multScen.maxRelRangeShift     = 3.5;             % maximum relative over and undershoot in % 
   multScen.rangeCombType        = 'combined';      % individual:  no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
   multScen.rangeGenType         = 'equidistant';   % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
   multScen.scenCombType         = 'individual';    % individual:  no combination of range and setup scenarios, 
                                                    % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                    % permuted:    create every possible combination of range and setup scenarios
   multScen.includeNomScen       = true;
   
% definition of scenarios for sampling
elseif ~pln.robOpt && pln.sampling
    requiredFields = {'numOfShiftScen','shiftSize','shiftGenType','shiftCombType',...
        'shiftGen1DIsotropy','numOfRangeShiftScen','maxAbsRangeShift','maxRelRangeShift',...
        'rangeCombType','rangeGenType','scenCombType','includeNomScen'};
  if sum(isfield(multScen,requiredFields)) == 0

     %% grid sampling
       % a) define shift scenarios
       multScen.numOfShiftScen       = [0 0 0];          % number of shifts in x y and z direction       
       multScen.shiftSize            = [3 3 3];          % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
       multScen.shiftGenType         = 'equidistant';    % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
       multScen.shiftCombType        = 'individual';     % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
                                                         % combined:    combine shift scenarios;                 number of shift scenarios is multScen.numOfShiftScen(1)
                                                         % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
       multScen.shiftGen1DIsotropy   = '+-';             % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

       % b) define range error scenarios                                                
       multScen.numOfRangeShiftScen  = 20; % number of absolute and/or relative range scnearios. 
                                                         % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
       multScen.maxAbsRangeShift     = 1;                % maximum absolute over and undershoot in mm   
       multScen.maxRelRangeShift     = 3.5;              % maximum relative over and undershoot in % 
       multScen.rangeCombType        = 'combined';       % individual: no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
       multScen.rangeGenType         = 'equidistant';    % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
       multScen.scenCombType         = 'individual';     % individual:  no combination of range and setup scenarios, 
                                                         % combined:    combine range and setup scenarios if their scenario number is consistent 
                                                         % permuted:    create every possible combination of range and setup scenarios
       multScen.includeNomScen       = false;
   
   %% random sampling
%    NumOfSamples = 20;
%    multScen.numOfShiftScen       = [NumOfSamples NumOfSamples NumOfSamples];  % number of shifts in x y and z direction       
%    multScen.shiftSize            = [3 3 3];          % maximum shift [mm]  % (e.g. prostate cases 5mm otherwise 3mm)
%    multScen.shiftGenType         = 'sampled';        % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
%    multScen.shiftCombType        = 'combined';       % individual:  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
%                                                      % combined:    combine shift scenarios;                 number of shift scenarios is multScen.numOfShiftScen(1)
%                                                      % permuted:    create every possible shift combination; number of shift scenarios is 8,27,64 ... 
%    multScen.shiftGen1DIsotropy   = '+-';             % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 
% 
%    % b) define range error scenarios                                                
%    multScen.numOfRangeShiftScen  = NumOfSamples;     % number of absolute and/or relative range scnearios. 
%                                                      % if absolute and relative range scenarios are defined then multScen.rangeCombType defines the resulting number of range scenarios
%    multScen.maxAbsRangeShift     = 1;                % maximum absolute over and undershoot in mm   
%    multScen.maxRelRangeShift     = 3.5;              % maximum relative over and undershoot in % 
%    multScen.rangeCombType        = 'combined';       % individual: no combination of absolute and relative range scenarios; combined:    combine absolute and relative range scenarios
%    multScen.rangeGenType         = 'sampled';        % equidistant: equidistant range shifts, sampled: sample range shifts from normal distribution
%    multScen.scenCombType         = 'combined';       % individual:  no combination of range and setup scenarios, 
%                                                      % combined:    combine range and setup scenarios if their scenario number is consistent 
%                                                      % permuted:    create every possible combination of range and setup scenarios
%    multScen.includeNomScen       = false;
%    
   
 elseif sum(isfield(multScen,requiredFields)) ~= numel(requiredFields)
       matRad_dispToConsole('Set of parameters is incomplete.',param,'error');
 end
   
else
   matRad_dispToConsole('matRad_setPlanUncertainties: Invalid combination',param,'error');
 end

%% create multiScen struct
pln.multScen           = matRad_setMultScen(multScen);
pln.numOfSamples       = pln.multScen.numOfScen;

%% get probabilities

pln.multScen.scenProb  = matRad_calcScenProb([0 0 0 0 0],[multScen.shiftSD multScen.rangeAbsSD multScen.rangeRelSD/100],...
                         pln.multScen.scenForProb,'probBins','normDist');
 
end


