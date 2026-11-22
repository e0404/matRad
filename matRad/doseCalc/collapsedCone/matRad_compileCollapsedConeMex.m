function success = matRad_compileCollapsedConeMex(dest, useOpenMP)
% matRad_compileCollapsedConeMex compiles the collapsed cone transport mex
% file for the current platform (MATLAB mex or Octave mkoctfile --mex)
%
% If the compilation with OpenMP fails, a serial version is compiled. With
% MATLAB on Linux the C++ runtime is linked statically to avoid conflicts
% with the runtime libraries shipped with MATLAB. Success is determined by
% calling the compiled mex file.
%
% call:
%   matRad_compileCollapsedConeMex()
%   matRad_compileCollapsedConeMex(dest)
%   matRad_compileCollapsedConeMex(dest, useOpenMP)
%   success = matRad_compileCollapsedConeMex(...)
%
% input:
%   dest:       (optional) destination folder of the mex file. Default:
%               location of this file
%   useOpenMP:  (optional) compile with OpenMP support. Default: true
%
% output:
%   success:    true if the compiled mex file can be called
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
if nargin < 2
    useOpenMP = true;
end

currDir = pwd;
cleaner = onCleanup(@() cd(currDir));
cd(dest);

success = false;
if useOpenMP
    matRad_compileCollapsedConeMexWithFlags(sourceFolder, true);
    success = matRad_testCollapsedConeMex();
    if ~success
        matRad_cfg.dispWarning('Compilation with OpenMP failed, compiling serial version.');
    end
end
if ~success
    matRad_compileCollapsedConeMexWithFlags(sourceFolder, false);
    success = matRad_testCollapsedConeMex();
end
if ~success
    matRad_cfg.dispWarning('Compilation of matRad_collapsedConeMex failed. The MATLAB reference implementation will be used.');
end

end

function matRad_compileCollapsedConeMexWithFlags(sourceFolder, useOpenMP)
matRad_cfg = MatRad_Config.instance();
sourceFile = fullfile(sourceFolder, 'matRad_collapsedConeMex.cpp');

if matRad_cfg.isOctave
    isMSVC = false;
else
    cc = mex.getCompilerConfigurations('C++', 'Selected');
    isMSVC = ~isempty(cc) && ~isempty(strfind(cc.ShortName, 'MSVC'));
end

if isMSVC
    compFlags = '/O2 /fp:precise';
    linkFlags = '';
    if useOpenMP
        compFlags = [compFlags ' /openmp'];
    end
else
    compFlags = '-std=c++11 -O2 -fPIC -ffp-contract=off';
    linkFlags = '';
    if useOpenMP
        compFlags = [compFlags ' -fopenmp'];
        linkFlags = '-fopenmp';
    end
    if matRad_cfg.isMatlab && isunix && ~ismac
        linkFlags = [linkFlags ' -static-libstdc++ -static-libgcc'];
    end
end

clear('matRad_collapsedConeMex');
try
    if matRad_cfg.isOctave
        oldFlags = {getenv('CXXFLAGS'), getenv('LDFLAGS')};
        restoreEnv = onCleanup(@() cellfun(@setenv, {'CXXFLAGS', 'LDFLAGS'}, oldFlags));
        setenv('CXXFLAGS', compFlags);
        setenv('LDFLAGS', linkFlags);
        mkoctfile('--mex', sourceFile);
    elseif isMSVC
        mex('-largeArrayDims', ['COMPFLAGS=$COMPFLAGS ' compFlags], sourceFile);
    else
        mex('-largeArrayDims', ['CXXFLAGS=$CXXFLAGS ' compFlags], ['LDFLAGS=$LDFLAGS ' linkFlags], sourceFile);
    end
catch err
    matRad_cfg.dispDebug('mex reported: %s\n', err.message);
end
end

function success = matRad_testCollapsedConeMex()
clear('matRad_collapsedConeMex');
try
    dose = matRad_collapsedConeMex(ones(2, 2, 2), ones(2, 2, 2), [1 1 1], [0; 0; 1], 1, 1, 0, 1);
    success = isequal(size(dose), [2 2 2]) && all(isfinite(dose(:)));
catch
    success = false;
end
end
