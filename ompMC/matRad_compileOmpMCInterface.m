function matRad_compileOmpMCInterface(dest,omcFolder)
% Compiles the ompMC interface (integrated as submodule)
%
% call
%   matRad_compileOmpMCInterface()
%   matRad_compileOmpMCInterface(dest)
%   matRad_compileOmpMCInterface(dest,sourceFolder)
%
% input:
%   dest:           (optional) destination for mex file. Default: location 
%                   of this file
%   sourceFolder:   (optional) path to ompMC . Default assumes its checked
%                   out in the submodules folder of matRad
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
matRad_cfg = MatRad_Config.instance();

env = matRad_getEnvironment();

if nargin < 1
    dest = fileparts(mfilename('fullpath'));
end

if nargin < 2
    omcFolder = [matRad_cfg.matRadRoot filesep 'submodules' filesep 'ompMC'];
end

sourceFolder = [omcFolder filesep 'src'];
interfaceFolder = [omcFolder filesep 'ucodes' filesep 'omc_matrad'];

mainFile = [interfaceFolder filesep 'omc_matrad.c'];

addFiles = {'ompmc.c','omc_utilities.c','omc_random.c'};
addFiles = cellfun(@(f) fullfile(sourceFolder,f),addFiles,'UniformOutput',false);

addFiles = strjoin(addFiles,' ');

if exist ('OCTAVE_VERSION','builtin')
    ccName = eval('mkoctfile -p CC');
else
    myCCompiler = mex.getCompilerConfigurations('C','Selected');
    ccName = myCCompiler.ShortName;
end

%These settings have only been tested for MSVC and g++. You may need to adapt for other compilers
if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
    flags{1,1} = 'COMPFLAGS';
    flags{1,2} = '/openmp';
    flags{2,1} = 'OPTIMFLAGS';
    flags{2,2} = '/O2';
else
    flags{1,1} = 'CFLAGS';
    flags{1,2} = '-std=gnu99 -fopenmp -O3';
    flags{2,1} = 'LDFLAGS';
    flags{2,2} = '-fopenmp';
    
end

includestring =  ['-I' sourceFolder];

flagstring = '';

%For Octave, the flags will be set in the environment, while they
%will be parsed as string arguments in MATLAB
for flag = 1:size(flags,1)
    if strcmp(env,'OCTAVE')
        preFlagContent = eval(['mkoctfile -p ' flags{flag,1}]);
        if ~isempty(preFlagContent)
            preFlagContent = preFlagContent(1:end-1); %Strip newline
        end
        newContent = [preFlagContent ' ' flags{flag,2}];
        setenv(flags{flag,1},newContent);
        matRad_cfg.dispDebug('Set compiler flag %s to %s\n',flags{flag,1},newContent);
    else
        flagstring = [flagstring flags{flag,1} '="$' flags{flag,1} ' ' flags{flag,2} '" '];
    end
end

mexCall = ['mex -largeArrayDims ' flagstring ' ' includestring ' ' mainFile ' ' addFiles];
matRad_cfg.dispDebug('Compiler call: %s\n',mexCall);

currDir = pwd;
cd(dest);
eval(mexCall);
cd(currDir);
end
