function options = matRad_probOptions(ct,cst,pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad APM options file
% 
% This script allows to define several options for analytical probabilistic
% modeling. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global INPUT_UNCERTAINTY;
global OUTPUT_UNCERTAINTY;
global SHOW_PARALLELTOOLBOX_WARNING;

SHOW_PARALLELTOOLBOX_WARNING = true;

options.probOpt.visBool       = 0;
options.probOpt.numOfWorkers  = 0; 

try
   [FlagParallToolBoxLicensed,msg]  = license('checkout','Distrib_Computing_Toolbox');
   if ~FlagParallToolBoxLicensed
      if SHOW_PARALLELTOOLBOX_WARNING
         matRad_dispToConsole(['Could not check out parallel computing toolbox. \n'],[],'warning');
         SHOW_PARALLELTOOLBOX_WARNING = false;
      end
   end
catch
   FlagParallToolBoxLicensed  = false;
end

% decide which version of variance calculation should run
options.probOpt.typeVarCalc    = 'parallel';   %parallel or serial 
options.probOpt.numOfWorkers   = 0;
    
if FlagParallToolBoxLicensed
    myCluster                    = parcluster('local');
    options.probOpt.numOfWorkers = myCluster.NumWorkers; 
    options.probOpt.typeVarCalc  = 'parallelParFor'; 
end

options.probOpt.omegaTarget  = false;               %  if true only target voxels are used for omega filling  

options.probOpt.LatCutOff    = 3.5^2;               % [sigma units] or set it to Inf if no additional cutoff should be used
                                                    % this option is only relevant for the std calculation

if strcmp(pln.radiationMode,'protons')
    options.probOpt.relDoseThreshold       = 0.01;  % std is only calculated for voxels having a greater relative dose than the threshold 
elseif strcmp(pln.radiationMode,'carbon')
    options.probOpt.relDoseThreshold       = 0.01;  % std is only calculated for voxels having a greater relative dose than the threshold 
end

options.probOpt.calcFullFrac        = false;       % calculate std for multiple number of fractions throughout the RT course
options.probOpt.calcFrac            = [1 5 10 30]; % calc std for the following fractionation scheme
options.probOpt.useReducedComp      = false;       % use reduced gaussian components 
options.probOpt.reducedComprelDose  = 0.02;        % if gauss component has a lower relative dose than the threshold then omitt it

% define which sources of uncertainty should be included
if ~isempty(INPUT_UNCERTAINTY)
   options.probOpt.InputUCT = INPUT_UNCERTAINTY;
else
   options.probOpt.InputUCT = 'phys'; %'phys', 'bio', 'biophys'
end

if ~isempty(OUTPUT_UNCERTAINTY)
   options.probOpt.OutputUCT = OUTPUT_UNCERTAINTY;
else
   options.probOpt.OutputUCT = 'phys'; %'phys', 'bio', 'biophys'
end

%%  define voxel list for which variance should be calculated
% use all voxels in the cst
%V = [cst{:,4}];
%options.probOpt.voxelList  = unique(vertcat(V{:}))'; 

% use iso center slice voxels
options.probOpt.slice        = round(pln.propStf.isoCenter(1,3)/ct.resolution.z);
options.probOpt.voxelList    = double(ct.cubeDim(1) * ct.cubeDim(2) * (options.probOpt.slice-1) + (1:1:(ct.cubeDim(1) * ct.cubeDim(2))));



end
