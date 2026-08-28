function matRad_concatenateRunfiles(varHelper, pathRunfiles)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function concatenates blocks for MCNP input to one runfile for each bixel
% calles 'MCNPrunfile_bixelN' where N is the bixel number.
% 
% call
%   matRad_concatenateRunfiles(varHelper, pathRunfiles)
%
% input
%   varHelper:      Helping variable with varHelper.totalNumberBixels and 
%                   varHelper.simPropMCNP.sourceBlockNames
%   pathRunfiles:   Path to previously generated MCNP runfiles
%
% output
%   none
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 12/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

oldPath = cd(pathRunfiles);
restoreDir = onCleanup(@() cd(oldPath));

if ismac || isunix
    for counterRunfile=1:varHelper.totalNumberBixels
        dummy_name = convertStringsToChars(varHelper.simPropMCNP.sourceBlockNames(counterRunfile));
        dummy_nameRunfile = ['MCNPrunfile_', num2str(counterRunfile),'bixel'];
        dummy_nameRunfile = convertStringsToChars(dummy_nameRunfile);
        system(['cat blockA.txt >> ', dummy_nameRunfile]);
        system(['cat blockB.txt >> ', dummy_nameRunfile]);
        system(['cat ', dummy_name,' >> ', dummy_nameRunfile]);
        system(['cat blockC_rest >> ', dummy_nameRunfile]);
    end
elseif ispc
    for counterRunfile=1:varHelper.totalNumberBixels
        dummy_name = convertStringsToChars(varHelper.simPropMCNP.sourceBlockNames(counterRunfile));
        dummy_nameRunfile = ['MCNPrunfile_', num2str(counterRunfile),'bixel'];
        dummy_nameRunfile = convertStringsToChars(dummy_nameRunfile);
        system(['type blockA.txt > ', dummy_nameRunfile]);
        system(['type blockB.txt >> ', dummy_nameRunfile]);
        system(['type ', dummy_name,' >> ', dummy_nameRunfile]);
        system(['type blockC_rest >> ', dummy_nameRunfile]);
    end
else
    matRad_cfg.dispInfo('Platform not supported but you can concatenate the blocks to one runfile by hand.\n');
end
