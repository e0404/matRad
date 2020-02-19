function matRad_compileOmpMCInterface()
% Compiles the ompMC interface (integrated as submodule)
%
% call
%   matRad_compileOmpMCInterface()
%
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

    if exist ('OCTAVE_VERSION','builtin')
      ccName = eval('mkoctfile -p CC');
    else
      myCCompiler = mex.getCompilerConfigurations('C','Selected');
      ccName = myCCompiler.ShortName;
    end


    %This needs to generalize better
    if ~isempty(strfind(ccName,'MSVC')) %Not use contains(...) because of octave
        flags{1,1} = 'COMPFLAGS';
        flags{1,2} = '$COMPFLAGS /openmp';
        flags{2,1} = 'OPTIMFLAGS';
        flags{2,2} = '$OPTIMFLAGS /O2';
        %flags = [optPrefix 'COMPFLAGS="$COMPFLAGS /openmp" ' optPrefix 'OPTIMFLAGS="$OPTIMFLAGS /O2"'];
    else
        flags{1,1} = 'CFLAGS';
        flags{1,2} = '--std=c99 -fopenmp -O2 -fPIC';
        flags{2,1} = 'LDFLAGS';
        flags{2,2} = '-fopenmp';
        %flags = [optPrefix 'CFLAGS="$CFLAGS -fopenmp -O2" ' optPrefix 'LDFLAGS="$LDFLAGS -fopenmp"'];
    end

    %flags = {};
    %flagstring = '-g ';
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
    
    mexCall = ['mex -largeArrayDims ' flagstring ' matRad_ompInterface.c'];
        
    disp(['Compiler call: ' mexCall]);
    eval(mexCall);
end
