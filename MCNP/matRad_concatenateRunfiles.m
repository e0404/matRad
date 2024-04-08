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

oldPath = pwd;
cd(pathRunfiles)

if ismac || isunix
    for counterRunfile=1:varHelper.totalNumberBixels
        dummy_name = convertStringsToChars(varHelper.simPropMCNP.sourceBlockNames(counterRunfile));
        dummy_nameRunfile = ['MCNPrunfile_bixel', num2str(counterRunfile)];
        dummy_nameRunfile = convertStringsToChars(dummy_nameRunfile);
        system(['cat blockA.txt >> ', dummy_nameRunfile]);
        system(['cat blockB.txt >> ', dummy_nameRunfile]);
        system(['cat ', dummy_name,' >> ', dummy_nameRunfile]);
        system(['cat blockC_rest >> ', dummy_nameRunfile]);
    end
elseif ispc
    for counterRunfile=1:varHelper.totalNumberBixels
        dummy_name = convertStringsToChars(varHelper.simPropMCNP.sourceBlockNames(counterRunfile));
        dummy_nameRunfile = ['MCNPrunfile_bixel', num2str(counterRunfile)];
        dummy_nameRunfile = convertStringsToChars(dummy_nameRunfile);
        system(['type blockA.txt > ', dummy_nameRunfile]);
        system(['type blockB.txt >> ', dummy_nameRunfile]);
        system(['type ', dummy_name,' >> ', dummy_nameRunfile]);
        system(['type blockC_rest >> ', dummy_nameRunfile]);
    end
else
    disp('Platform not supported but you can concatenate the blocks to one runfile by hand.')
end

cd(oldPath);