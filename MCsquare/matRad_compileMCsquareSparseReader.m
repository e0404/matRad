function matRad_compileMCsquareSparseReader()
% Compiles the sparse mcsquare reader as mex interface for the current
% platform
%
% call
%   matRad_compileMCsquareSparseReader()
%
%
% References
%
%   https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.4943377
%   http://www.openmcsquare.org/
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

if exist ('OCTAVE_VERSION', 'builtin')
    ccName = eval('mkoctfile -p CXX');
else
    myCCompiler = mex.getCompilerConfigurations('C','Selected');
    ccName = myCCompiler.ShortName;
end

%This needs to generalize better
if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
    flags{1,1} = 'COMPFLAGS';
    flags{1,2} = '/O2';
else
    flags{1,1} = 'CXXFLAGS';
    flags{1,2} = '-std=c++11 -O2 -fPIC';
end

%flags = {};
%flagstring = '-g ';
flagstring = '';

%For Octave, the flags will be set in the environment, while they
%will be parsed as string arguments in MATLAB
for flag = 1:size(flags,1)
    if exist ('OCTAVE_VERSION', 'builtin')
        setenv(flags{flag,1},flags{flag,2});
    else
        flagstring = [flagstring flags{flag,1} '="' flags{flag,2} '" '];
    end
end

mexCall = ['mex -largeArrayDims ' flagstring ' matRad_sparseBeamletsReaderMCsquare.cpp'];

disp(['Compiler call: ' mexCall]);
eval(mexCall);


end
