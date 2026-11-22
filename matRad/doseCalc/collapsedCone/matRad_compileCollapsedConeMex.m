function success = matRad_compileCollapsedConeMex(dest, useOpenMP, staticRuntime)
% matRad_compileCollapsedConeMex compiles the collapsed cone transport mex
% file for the current platform (MATLAB mex or Octave mkoctfile --mex)
%
% After compilation the mex file is loaded and called on a small problem.
% Compilation errors and mex files that cannot be called raise an error.
% The only tolerated compiler message is a spurious "is not a MEX file"
% report issued by some MATLAB releases for binaries that load and run.
%
% call:
%   matRad_compileCollapsedConeMex()
%   matRad_compileCollapsedConeMex(dest)
%   matRad_compileCollapsedConeMex(dest, useOpenMP)
%   matRad_compileCollapsedConeMex(dest, useOpenMP, staticRuntime)
%   success = matRad_compileCollapsedConeMex(...)
%
% input:
%   dest:           (optional) destination folder of the mex file.
%                   Default: location of this file
%   useOpenMP:      (optional) compile with OpenMP support. Default: true
%   staticRuntime:  (optional) link libstdc++ and libgcc statically
%                   (MATLAB on Linux only). Needed when the system compiler
%                   is newer than the one supported by the running MATLAB
%                   release and the mex file conflicts with the C++ runtime
%                   shipped with MATLAB. Default: false
%
% output:
%   success:        true if the compiled mex file can be called
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
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

sourceFolder = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(dest)
    dest = sourceFolder;
end
if nargin < 2 || isempty(useOpenMP)
    useOpenMP = true;
end
if nargin < 3 || isempty(staticRuntime)
    staticRuntime = false;
end

if staticRuntime && ~(matRad_cfg.isMatlab && isunix && ~ismac)
    matRad_cfg.dispWarning('Static C++ runtime linking is only supported for MATLAB on Linux and will be ignored.');
    staticRuntime = false;
end

[compiler, flags] = matRad_collapsedConeMexFlags(useOpenMP, staticRuntime);
sourceFile = fullfile(sourceFolder, 'matRad_collapsedConeMex.cpp');

currDir = pwd;
cleaner = onCleanup(@() cd(currDir));
cd(dest);
clear('matRad_collapsedConeMex');

try
    switch compiler
        case 'octave'
            oldFlags = {getenv('CXXFLAGS'), getenv('LDFLAGS')};
            restoreEnv = onCleanup(@() cellfun(@setenv, {'CXXFLAGS', 'LDFLAGS'}, oldFlags));
            setenv('CXXFLAGS', flags.compile);
            setenv('LDFLAGS', flags.link);
            mkoctfile('--mex', sourceFile);
        case 'msvc'
            mex('-largeArrayDims', ['COMPFLAGS=$COMPFLAGS ' flags.compile], sourceFile);
        otherwise
            mex('-largeArrayDims', ['CXXFLAGS=$CXXFLAGS ' flags.compile], ['LDFLAGS=$LDFLAGS ' flags.link], sourceFile);
    end
catch err
    if isempty(strfind(err.message, 'is not a MEX file'))
        matRad_cfg.dispError('Compilation of matRad_collapsedConeMex failed:\n%s', err.message);
    end
    matRad_cfg.dispDebug('mex reported "%s", checking the compiled binary.\n', err.message);
end

[success, msg] = matRad_testCollapsedConeMex();
if ~success
    matRad_cfg.dispError(['The compiled matRad_collapsedConeMex cannot be called: %s\n' ...
                          'If your compiler is newer than the one supported by MATLAB, try staticRuntime = true.'], msg);
end

end

function [compiler, flags] = matRad_collapsedConeMexFlags(useOpenMP, staticRuntime)
matRad_cfg = MatRad_Config.instance();

if matRad_cfg.isOctave
    compiler = 'octave';
else
    cc = mex.getCompilerConfigurations('C++', 'Selected');
    if isempty(cc)
        matRad_cfg.dispError('No C++ compiler configured for mex. Run ''mex -setup C++''.');
    end
    if ~isempty(strfind(cc.ShortName, 'MSVC'))
        compiler = 'msvc';
    else
        compiler = 'gcc';
    end
end

if strcmp(compiler, 'msvc')
    flags.compile = '/O2 /fp:precise';
    flags.link = '';
    if useOpenMP
        flags.compile = [flags.compile ' /openmp'];
    end
    return
end

flags.compile = '-std=c++11 -O2 -fPIC -ffp-contract=off';
flags.link = '';
if useOpenMP
    flags.compile = [flags.compile ' -fopenmp'];
    flags.link = '-fopenmp';
end
if staticRuntime
    flags.link = strtrim([flags.link ' -static-libstdc++ -static-libgcc']);
end
end

function [success, msg] = matRad_testCollapsedConeMex()
clear('matRad_collapsedConeMex');
msg = '';
try
    dose = matRad_collapsedConeMex(ones(2, 2, 2), ones(2, 2, 2), [1 1 1], [0; 0; 1], 1, 1);
    success = isequal(size(dose), [2 2 2]) && all(isfinite(dose(:)));
    if ~success
        msg = 'unexpected result of the test call';
    end
catch err
    success = false;
    msg = err.message;
end
end
