function matRad_compileOmpMCInterface(dest,omcFolder)
% Compiles the ompMC interface (integrated as submodule)
%
% call
%   matRad_compileOmpMCInterface()
%   matRad_compileOmpMCInterface(dest)
%   matRad_compileOmpMCInterface(dest,sourceFolder)
%
% input:
%   dest:           destination for mex file. Default: location of this
%                   file
%   sourceFolder:   path to ompMC (optional). Default assumes its checked
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

global matRad_cfg; 
matRad_cfg = MatRad_Config.instance();

if nargin < 1
    dest = fileparts(mfilename('fullpath'));
end

if nargin < 2
    omcFolder = [matRad_cfg.matRadRoot filesep 'submodules' filesep 'ompMC'];
end

sourceFolder = [omcFolder filesep 'src'];
interfaceFolder = [omcFolder filesep 'ucodes' filesep 'omc_matrad'];

mainFile = [interfaceFolder filesep 'omc_matRad.c'];

addFiles = {'ompmc.c','omc_utilities.c','omc_random.c'};
addFiles = cellfun(@(f) fullfile(sourceFolder,f),addFiles,'UniformOutput',false);

addFiles = strjoin(addFiles,' ');

if exist ('OCTAVE_VERSION','builtin')
    ccName = eval('mkoctfile -p CC');
else
    myCCompiler = mex.getCompilerConfigurations('C','Selected');
    ccName = myCCompiler.ShortName;
end


%This needs to generalize better
if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
    flags{1,1} = 'COMPFLAGS';
    flags{1,2} = ['/I' sourceFolder ' $COMPFLAGS /openmp'];
    flags{2,1} = 'OPTIMFLAGS';
    flags{2,2} = '$OPTIMFLAGS /O2';
else
    flags{1,1} = 'CFLAGS';
    flags{1,2} = ['-I' sourceFolder ' -std=c99 -fopenmp -O2 -fPIC'];
    flags{2,1} = 'LDFLAGS';
    flags{2,2} = '-fopenmp';
end

flagstring = '';

%For Octave, the flags will be set in the environment, while they
%will be parsed as string arguments in MATLAB
for flag = 1:size(flags,1)
    if exist ('OCTAVE_VERSION','builtin')
        setenv(flags{flag,1},flags{flag,2});
    else
        flagstring = [flagstring flags{flag,1} '="' flags{flag,2} '" '];
    end
end

mexCall = ['mex -largeArrayDims ' flagstring ' ' mainFile ' ' addFiles];
matRad_cfg.dispDebug('Compiler call: %s\n',mexCall);

currDir = pwd;
cd(dest);
eval(mexCall);
cd(currDir);
end
