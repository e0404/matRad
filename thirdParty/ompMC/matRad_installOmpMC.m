function installDir = matRad_installOmpMC(varargin)
% matRad_installOmpMC installs the ompMC photon Monte Carlo engine
%
% ompMC is licensed under the GNU General Public License v3, matRad under the
% 3-clause BSD license. matRad therefore does not ship ompMC - neither its
% sources nor a compiled mex interface - and instead obtains it on demand, on
% your machine, after you accepted the GPL. This function is that step.
%
% It either downloads the pre-compiled mex interface of a pinned ompMC release
% matching your platform, or builds one from the ompMC sources in the
% submodules folder. Everything ends up in thirdParty/ompMC, which is ignored
% by git.
%
% Whichever way it got there, the interface is called once before the install
% counts as done. A pre-compiled release that does not load - built against a
% newer MATLAB than yours, say - is therefore not left behind as a broken
% installation: it is replaced by one built here.
%
% call:
%   matRad_installOmpMC()
%   matRad_installOmpMC('acceptLicense',true)
%   matRad_installOmpMC('mode','build')
%   matRad_installOmpMC('uninstall',true)
%
% input (name-value pairs):
%   acceptLicense:  accept the GPL of ompMC without being asked. Required for
%                   non-interactive use. Default: false (you will be asked)
%   mode:           'auto'     download a release if one exists for this
%                              platform, and build from source if there is
%                              none or if what was downloaded turns out not
%                              to load here (default)
%                   'download' download only, fail if there is no release or
%                              if it does not load
%                   'build'    build from source only
%   sourceFolder:   ompMC sources to build from.
%                   Default: <matRadRoot>/submodules/ompMC
%   force:          re-install even if ompMC is already installed.
%                   Default: false
%   uninstall:      remove an installed ompMC instead of installing it.
%                   Default: false
%
% output:
%   installDir:     folder ompMC was installed to (empty when uninstalling)
%
% References
%   https://github.com/e0404/ompMC
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

release = matRad_ompMCrelease();

p = inputParser;
p.addParameter('acceptLicense', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('mode', 'auto', @(x) any(strcmpi(x, {'auto', 'download', 'build'})));
p.addParameter('sourceFolder', fullfile(matRad_cfg.matRadRoot, 'submodules', 'ompMC'), @ischar);
p.addParameter('force', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('uninstall', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

installDir = fullfile(matRad_cfg.matRadRoot, 'thirdParty', 'ompMC');
mode = lower(p.Results.mode);

%% Uninstall
if p.Results.uninstall
    matRad_ompMCclean(installDir);
    matRad_cfg.dispInfo('ompMC removed from %s.\n', installDir);
    installDir = '';
    return
end

%% Nothing to do?
% Asked of the engine rather than answered here, so that whatever it refuses
% to load - a binary for another octave, say - counts as needing installing
if DoseEngines.matRad_PhotonOmpMCEngine.checkInstallation() && ~p.Results.force
    matRad_cfg.dispInfo('ompMC is already installed in %s. Pass ''force'',true to re-install.\n', installDir);
    return
end

% A mex file that has been called once cannot be overwritten while it is
% loaded - ompMC locks itself in memory because an OpenMP mex file cannot
% safely be unloaded once a parallel region has run.
if matRad_ompMCisLoaded()
    matRad_cfg.dispError(['ompMC has already been loaded in this session and cannot be replaced. ' ...
                          'Please restart %s and run matRad_installOmpMC again.'], matRad_getEnvironment());
end

%% License
if ~all(p.Results.acceptLicense)
    matRad_ompMCaskForLicense(release);
end

%% Install
asset = matRad_ompMCassetForThisPlatform(release);

switch mode
    case 'download'
        if isempty(asset)
            matRad_cfg.dispError(['ompMC %s has no pre-compiled release for this platform. ' ...
                                  'Use matRad_installOmpMC(''mode'',''build'') to build it from source.'], release.tag);
        end
        doDownload = true;
    case 'build'
        doDownload = false;
    otherwise
        doDownload = ~isempty(asset);
        if ~doDownload
            matRad_cfg.dispInfo('No pre-compiled ompMC release for this platform, building from source instead.\n');
        end
end

if doDownload
    matRad_ompMCclean(installDir);
    matRad_ompMCdownload(release, asset, installDir);

    loadError = matRad_ompMCfinalize(installDir);

    if isempty(loadError)
        return
    end

    % A package that does not load is not a package for this machine,
    % whatever it was built for - too new a MATLAB is one way, and the
    % message it fails with rarely says which. Building says it in terms of
    % the compiler and the MATLAB that are actually here.
    if strcmp(mode, 'download')
        matRad_cfg.dispError(['The pre-compiled ompMC %s could not be loaded: %s\n' ...
                              'Build it for this machine with matRad_installOmpMC(''mode'',''build'').'], ...
                             asset.name, loadError);
    end

    matRad_cfg.dispWarning(['The pre-compiled ompMC %s could not be loaded (%s), ' ...
                            'building it from source instead.'], asset.name, loadError);
end

matRad_ompMCclean(installDir);
matRad_ompMCbuild(release, p.Results.sourceFolder, installDir);

loadError = matRad_ompMCfinalize(installDir);

if ~isempty(loadError)
    matRad_cfg.dispError('ompMC was built in %s, but the mex interface could not be called: %s', ...
                         installDir, loadError);
end

end

%% ------------------------------------------------------------------------
function loadError = matRad_ompMCfinalize(installDir)
% Completes an installation and reports whether what it put there can
% actually be called, as the message it failed with or empty if it can.

matRad_cfg = MatRad_Config.instance();

% So that a later session can tell whether it is able to load what is here -
% a mex file for another octave has the same name and fails only when called
DoseEngines.matRad_PhotonOmpMCEngine.writeBinaryTag();

% ompMC writes its own diagnostics here, and the release packages do not
% carry the (empty) folder
outputDir = fullfile(installDir, 'output');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% The mex file sits in bin/, which is only on the path from the next
% matRad_rc on, so this session gets it now
addpath(fullfile(installDir, 'bin'));

% A mex file that failed to load once is remembered as broken, which would
% hide a working one put in its place afterwards
clear('omc_matrad');
rehash;

loadError = '';

try
    installedVersion = omc_matrad('version');
catch ME
    loadError = ME.message;
    return
end

matRad_cfg.dispInfo('ompMC %s installed in %s\n', installedVersion, installDir);
end

%% ------------------------------------------------------------------------
function release = matRad_ompMCrelease()
% The ompMC release matRad is built against. Everything that has to change
% when it is bumped lives here: the tag, and the SHA-256 of every asset as
% published with the release.

release.tag     = 'v0.3.0';
release.version = '0.3.0';
release.repoUrl = 'https://github.com/e0404/ompMC';

% platform (as in the asset name), SHA-256 of ompmc-<version>-<platform>.zip
release.assets = { ...
                  'linux-arm64',             '609933e4a7ac9426e8afe7165e13eecc9a5c0d9d21a05a5662d81a5279760455'; ...
                  'linux-x64',               '3fabb7b059373a5125c92afc09e2cdd7400a46323e85884deb30526e9c1571db'; ...
                  'linux-x64-octave8',       '068090d7e9d382984495cc65a56645b216864ac9a0148ca5f21a96623b592aa8'; ...
                  'macos-arm64',             'b5bcfc03e16b06bd4e07527d306dd4693baea5c2d262fcbbd9ee716d3efaf171'; ...
                  'macos-x64',               'd86b391c4df6480cb085a5662d50dd1ac54ce57238c3870853cf228e4e90ba07'; ...
                  'windows-x64-llvm-mingw',  '6bf446e6e0dc778a3e7ebd34ae15adc324113a50f070eb9892693b63c1ef8bfe'; ...
                  'windows-x64-mingw',       '80c7be21c9702a9bc797b2dc8f3d3cd35fc9b510e8bd8ceb0c021ab0b756e861'; ...
                  'windows-x64-msvc',        'ccb6cbc2cd607a3f3fa87afd9058a70c8205c067ae76aebc4f186ca01e07567c'; ...
                  'windows-x64-octave8',     '426ab6d9c3d065f5e3c7d38ebba702e054e6aa4a8e3d7d4fb7325b5480cff4e4'; ...
                  'windows-x64-octave10',    '93a27f3542c3f7cb8f19bffb30467da884ffcad41be2d1f2fab8c812b6f38707'};
end

%% ------------------------------------------------------------------------
function asset = matRad_ompMCassetForThisPlatform(release)
% Picks the release asset for the running environment, or returns empty if
% the release does not have one.

isArm = matRad_ompMCisArm();

% ompMC packages one asset per ABI bucket, which is what the engine names to
% decide whether it can load a given binary - the same question, asked before
% rather than after the download
platform = DoseEngines.matRad_PhotonOmpMCEngine.binaryTag();

if ~strcmp(platform, 'matlab')
    if strcmp(platform, 'octave')
        % an octave whose version did not parse
        asset = '';
        return
    end

    if ispc
        platform = ['windows-x64-' platform];
    elseif ismac || isArm
        % no macOS and no arm64 Octave packages so far
        platform = '';
    else
        platform = ['linux-x64-' platform];
    end
else
    if ispc
        % The MinGW builds carry their own runtime dlls, the MSVC one only
        % needs the Visual C++ redistributable, which any machine running
        % MATLAB has
        platform = 'windows-x64-msvc';
    elseif ismac && isArm
        platform = 'macos-arm64';
    elseif ismac
        platform = 'macos-x64';
    elseif isArm
        platform = 'linux-arm64';
    else
        platform = 'linux-x64';
    end
end

asset = '';

if isempty(platform)
    return
end

ix = find(strcmp(release.assets(:, 1), platform), 1);

if ~isempty(ix)
    asset.platform = platform;
    asset.name     = sprintf('ompmc-%s-%s.zip', release.version, platform);
    asset.sha256   = release.assets{ix, 2};
    asset.url      = sprintf('%s/releases/download/%s/%s', release.repoUrl, release.tag, asset.name);
end
end

%% ------------------------------------------------------------------------
function isArm = matRad_ompMCisArm()
% Both environments spell the architecture differently, and neither spells
% it the same way on every system

arch = lower(computer());

isArm = ~isempty(strfind(arch, 'arm')) || ~isempty(strfind(arch, 'aarch64')) || strcmp(arch, 'maca64'); %#ok<STREMP>
end

%% ------------------------------------------------------------------------
function matRad_ompMCaskForLicense(release)
% Asks whether to install ompMC, which is also where its GPL is accepted -
% one question rather than two, since the answer is the same one. Skipped by
% passing 'acceptLicense', which is how a script says it in advance.

matRad_cfg = MatRad_Config.instance();

notice = { ...
          sprintf('ompMC %s is not part of matRad and is not distributed with it.', release.tag), ...
          '', ...
          'It is free software by Edgardo Doerner and contributors, licensed under', ...
          'the GNU General Public License v3 - not under the 3-clause BSD license', ...
          'of matRad. Installing it combines the two on this machine, which is your', ...
          'right to do; passing the result on to anyone else means passing on the', ...
          'terms of the GPL with it.', ...
          '', ...
          sprintf('    %s', release.repoUrl), ...
          sprintf('    %s/blob/%s/LICENSE', release.repoUrl, release.tag), ...
          '', ...
          'Install ompMC now, accepting the GNU General Public License v3 for it?'};

if matRad_ompMCcanShowDialog()
    answer = questdlg(notice, sprintf('Install ompMC %s?', release.tag), ...
                      'Accept and install', 'Cancel', 'Cancel');
    accepted = strcmp(answer, 'Accept and install');
else
    fprintf('\n  %s\n\n', strjoin(notice, sprintf('\n  ')));

    try
        answer = input('  [y/N] ', 's');
    catch
        % a session that cannot be asked at all, e.g. matlab -batch
        answer = '';
    end

    fprintf('\n');
    accepted = any(strcmpi(strtrim(answer), {'y', 'yes'}));
end

if accepted
    return
end

% Somebody who said no knows why this failed. An empty answer is the other
% case: a session that could not be asked, which is worth telling how to
% answer in advance.
if isempty(strtrim(answer))
    matRad_cfg.dispError(['ompMC was not installed because its license was not accepted. ' ...
                          'This session could not be asked - accept it beforehand with ' ...
                          'matRad_installOmpMC(''acceptLicense'',true).']);
else
    matRad_cfg.dispError('ompMC was not installed because its license was not accepted.');
end
end

%% ------------------------------------------------------------------------
function tf = matRad_ompMCcanShowDialog()
% Whether this session can put a dialog in front of somebody. A session that
% cannot is not necessarily one that cannot be asked at all - a console still
% takes an answer - so this only decides how to ask, not whether to.

matRad_cfg = MatRad_Config.instance();

tf = false;

if matRad_cfg.disableGUI
    return
end

if matRad_cfg.isOctave
    % true only with the octave GUI, not for octave-cli
    tf = exist('isguirunning', 'builtin') == 5 && isguirunning();
else
    % false under -batch and -nodisplay, where a dialog would either fail or
    % wait for a click that cannot come
    tf = usejava('desktop');
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCdownload(release, asset, installDir)
% Downloads and unpacks a pre-compiled ompMC release

matRad_cfg = MatRad_Config.instance();

% Unpacked next to where it is going rather than in the system temp folder,
% so that moving it into place stays within one volume - see
% matRad_ompMCmoveContent. Nothing below installDir is tracked, and it is
% removed either way.
tmpDir = fullfile(installDir, 'download');
matRad_ompMCrmdir(tmpDir);
mkdir(tmpDir);

cleanTmp = onCleanup(@() matRad_ompMCrmdir(tmpDir)); %#ok<NASGU>

zipFile = fullfile(tmpDir, asset.name);

matRad_cfg.dispInfo('Downloading %s...\n', asset.url);

try
    matRad_ompMCfetch(asset.url, zipFile);
catch ME
    matRad_cfg.dispError('Could not download ompMC from %s:\n%s', asset.url, ME.message);
end

matRad_ompMCverifyChecksum(zipFile, asset.sha256);

matRad_cfg.dispInfo('Unpacking...\n');
matRad_ompMCunzip(zipFile, tmpDir);

% The archive holds a single folder named like itself, whose content is what
% we want in the install folder
[~, assetBase] = fileparts(asset.name);
contentDir = fullfile(tmpDir, assetBase);

if ~exist(contentDir, 'dir')
    matRad_cfg.dispError('Unexpected layout in %s: no folder %s.', asset.name, assetBase);
end

matRad_ompMCmoveContent(contentDir, installDir);

% The release packages do not carry the license text, so it is fetched from
% the tag the binaries were built from
licenseUrl = sprintf('%s/raw/%s/LICENSE', release.repoUrl, release.tag);

try
    matRad_ompMCfetch(licenseUrl, fullfile(installDir, 'LICENSE'));
catch
    matRad_cfg.dispWarning('Could not download the ompMC license text from %s.', licenseUrl);
end

matRad_ompMCwriteOrigin(installDir, release, sprintf('pre-compiled release asset %s (sha256 %s)', asset.name, asset.sha256));
end

%% ------------------------------------------------------------------------
function matRad_ompMCunzip(zipFile, destDir)
% Unpacks a zip archive. Octave's unzip hands the work to an unzip program
% that a Windows machine does not usually have, so where the built-in one
% cannot do it, the system is asked directly - through bsdtar, which reads
% zip files and comes with Windows, or through unzip everywhere else.

try
    unzip(zipFile, destDir);
    return
catch ME
    builtinError = ME.message;
end

if ispc
    unpackCall = sprintf('tar -xf "%s" -C "%s"', zipFile, destDir);
else
    unpackCall = sprintf('unzip -o -q "%s" -d "%s"', zipFile, destDir);
end

[status, msg] = system(unpackCall);

if status ~= 0
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Could not unpack %s.\nBuilt-in unzip said: %s\n%s said: %s', ...
                         zipFile, builtinError, unpackCall, msg);
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCfetch(url, target)
% Downloads a file with whatever the environment has. MATLAB has websave;
% octave has urlwrite, but from octave 11 on that gives up after five
% seconds, which a release package does not download in, so where weboptions
% can lift the limit webread is asked instead. It hands the bytes back as
% char, one per byte - and if that ever stopped being true, the checksum the
% caller verifies afterwards is what would say so.

if exist('websave', 'file') == 2 || exist('websave', 'builtin') == 5
    websave(target, url);
    return
end

if exist('weboptions', 'file') == 2 || exist('weboptions', 'builtin') == 5
    data = webread(url, weboptions('Timeout', 300));

    fid = fopen(target, 'wb');

    if fid < 0
        error('matRad:ompMC:download', 'Could not write %s.', target);
    end

    fwrite(fid, data, 'uint8');
    fclose(fid);
    return
end

[~, success, msg] = urlwrite(url, target);

if ~success
    error('matRad:ompMC:download', '%s', msg);
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCverifyChecksum(fileName, expected)

matRad_cfg = MatRad_Config.instance();

actual = matRad_ompMCsha256(fileName);

if isempty(actual)
    matRad_cfg.dispWarning(['Could not compute a SHA-256 checksum in this environment, ' ...
                            'the downloaded ompMC release was not verified.']);
    return
end

if ~strcmpi(actual, expected)
    matRad_cfg.dispError(['Checksum mismatch for %s.\nexpected: %s\nfound:    %s\n' ...
                          'The download is corrupt or does not come from the pinned ompMC release, ' ...
                          'and was not installed.'], fileName, expected, actual);
end

matRad_cfg.dispDebug('Checksum of %s verified.\n', fileName);
end

%% ------------------------------------------------------------------------
function sha = matRad_ompMCsha256(fileName)
% SHA-256 of a file. Returns empty if neither environment offers one, which
% only happens for MATLAB started without a JVM.

sha = '';

fid = fopen(fileName, 'rb');

if fid < 0
    return
end

bytes = fread(fid, inf, '*uint8');
fclose(fid);

if strcmp(matRad_getEnvironment(), 'OCTAVE')
    sha = hash('SHA256', char(bytes'));
    return
end

if ~usejava('jvm')
    return
end

digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);

% digest() gives signed bytes, which have to be read as unsigned before they
% can be printed as the usual hex string
sha = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

%% ------------------------------------------------------------------------
function matRad_ompMCbuild(release, sourceFolder, installDir)
% Builds the ompMC mex interface from source, preferring CMake - which is how
% ompMC builds itself - and falling back to a plain mex call.

matRad_cfg = MatRad_Config.instance();

if ~exist(fullfile(sourceFolder, 'CMakeLists.txt'), 'file')
    matRad_cfg.dispError(['No ompMC sources in %s. Check out the submodule with\n' ...
                          '    git submodule update --init submodules/ompMC\n' ...
                          'or point matRad_installOmpMC at a source folder with ''sourceFolder''.'], sourceFolder);
end

mexFile = '';

% An octave mex file has to be built with the toolchain octave itself was
% built with, which is precisely what mkoctfile is. On Windows that is MinGW,
% while CMake picks up MSVC and ompMC then refuses to build the octave target
% at all - so there it is mkoctfile's job from the start.
useCMake = matRad_ompMChasCMake() && ~(matRad_cfg.isOctave && ispc);

if useCMake
    try
        mexFile = matRad_ompMCbuildWithCMake(sourceFolder);
    catch ME
        matRad_cfg.dispWarning('Building ompMC with CMake failed (%s), falling back to a plain mex call.', ME.message);
    end
else
    matRad_cfg.dispInfo('Building ompMC with a plain mex call.\n');
end

if isempty(mexFile)
    mexFile = matRad_ompMCbuildWithMex(sourceFolder);
end

% Assemble the same layout a release package has
binDir = fullfile(installDir, 'bin');
mkdir(binDir);
copyfile(mexFile, binDir);

for dataFolder = {'data', 'pegs4', 'spectra'}
    src = fullfile(sourceFolder, dataFolder{1});

    if ~exist(src, 'dir')
        matRad_cfg.dispError('The ompMC sources in %s do not contain the %s folder.', sourceFolder, dataFolder{1});
    end

    matRad_ompMCcopyDir(src, fullfile(installDir, dataFolder{1}));
end

copyfile(fullfile(sourceFolder, 'LICENSE'), installDir);

matRad_ompMCwriteOrigin(installDir, release, sprintf('built from the sources in %s', sourceFolder));
end

%% ------------------------------------------------------------------------
function matRad_ompMCcopyDir(fromDir, toDir)
% Copies a directory's content, one entry at a time into a destination that
% exists. Copying a directory onto a name that does not exist yet leaves
% octave's copyfile asking Windows whether that name means a file, and a
% wildcard is passed on literally rather than expanded.

matRad_cfg = MatRad_Config.instance();

if ~exist(toDir, 'dir')
    mkdir(toDir);
end

entries = dir(fromDir);

for i = 1:numel(entries)
    if any(strcmp(entries(i).name, {'.', '..'}))
        continue
    end

    from = fullfile(fromDir, entries(i).name);

    if entries(i).isdir
        matRad_ompMCcopyDir(from, fullfile(toDir, entries(i).name));
        continue
    end

    [success, msg] = copyfile(from, toDir);

    if ~success
        matRad_cfg.dispError('Could not copy %s into %s: %s', from, toDir, msg);
    end
end
end

%% ------------------------------------------------------------------------
function tf = matRad_ompMChasCMake()

[status, ~] = system('cmake --version');
tf = (status == 0);
end

%% ------------------------------------------------------------------------
function mexFile = matRad_ompMCbuildWithCMake(sourceFolder)
% Builds through ompMC's own CMake project, which is the configuration its
% releases are built in

matRad_cfg = MatRad_Config.instance();

buildDir = fullfile(tempname, 'ompMC-build');
mkdir(buildDir);

configureCall = sprintf('cmake -S "%s" -B "%s" -DCMAKE_BUILD_TYPE=Release -DOMPMC_BUILD_DOSXYZ=OFF', sourceFolder, buildDir);

if strcmp(matRad_getEnvironment(), 'OCTAVE')
    % Point it at the octave that is running, the way MATLAB is pointed at
    % below - CMake finds octave through octave-config on the PATH, which an
    % octave started by its full path is not on
    configureCall = sprintf('%s -DOMPMC_BUILD_MATRAD_MEX=OFF -DOMPMC_BUILD_MATRAD_OCT=ON -DOctave_ROOT="%s"', ...
                            configureCall, OCTAVE_HOME());
else
    % Build against the MATLAB that is running, not whichever one CMake finds
    configureCall = sprintf('%s -DOMPMC_BUILD_MATRAD_MEX=ON -DOMPMC_BUILD_MATRAD_OCT=OFF -DMatlab_ROOT_DIR="%s"', configureCall, matlabroot);
end

matRad_cfg.dispInfo('Configuring ompMC...\n');
matRad_cfg.dispDebug('%s\n', configureCall);
[status, msg] = system(configureCall);

if status ~= 0
    error('cmake configure failed:\n%s', msg);
end

buildCall = sprintf('cmake --build "%s" --config Release --parallel', buildDir);

matRad_cfg.dispInfo('Building ompMC...\n');
matRad_cfg.dispDebug('%s\n', buildCall);
[status, msg] = system(buildCall);

if status ~= 0
    error('cmake build failed:\n%s', msg);
end

mexFile = matRad_ompMCfindMexFile(fullfile(buildDir, 'bin'));

if isempty(mexFile)
    error('the build produced no omc_matrad mex file in %s', fullfile(buildDir, 'bin'));
end
end

%% ------------------------------------------------------------------------
function mexFile = matRad_ompMCbuildWithMex(sourceFolder)
% Fallback for machines without CMake: hand the whole core library and the
% matRad user code to mex in one call.
%
% ompMC's own build generates omc_version.h from the version in its
% CMakeLists, so that is done here as well. What CMake compiles conditionally
% - the Octave SOVERSION marker, which mkoctfile emits by itself, the MSVC
% export definition, which mex handles, and the Intel OpenMP shim - is left
% out.

matRad_cfg = MatRad_Config.instance();

env = matRad_getEnvironment();

srcDir       = fullfile(sourceFolder, 'src');
ucodeDir     = fullfile(sourceFolder, 'ucodes', 'omc_matrad');
generatedDir = fullfile(tempname, 'ompMC-generated');

mkdir(generatedDir);

matRad_ompMCwriteVersionHeader(sourceFolder, generatedDir);

% Everything in src/ is what ompMC's own build makes its core library from,
% and src/compat is not globbed because only the command line user code needs
% it. mex is called in function form, so none of these need quoting.
sourceFiles = dir(fullfile(srcDir, '*.c'));
sourceFiles = cellfun(@(f) fullfile(srcDir, f), {sourceFiles.name}, 'UniformOutput', false);
sourceFiles = [{fullfile(ucodeDir, 'omc_matrad.c')} sourceFiles];

% These settings have only been tested for MSVC and gcc. You may need to
% adapt them for other compilers
if strcmp(env, 'OCTAVE')
    ccName = evalc('mkoctfile -p CC');
else
    myCCompiler = mex.getCompilerConfigurations('C', 'Selected');
    ccName = myCCompiler.ShortName;
end

if ~isempty(strfind(ccName, 'MSVC')) %#ok<STREMP> - contains() does not exist in Octave
    flags = {'COMPFLAGS', '/openmp'; 'OPTIMFLAGS', '/O2'};
elseif strcmp(env, 'OCTAVE')
    % Only CFLAGS, which mkoctfile also passes when it links. Its LDFLAGS
    % carry octave's own library paths, and appending to those means handing
    % them back unquoted - which breaks every octave installed somewhere like
    % C:\Program Files.
    flags = {'CFLAGS', '-std=gnu99 -fopenmp -O3'};
else
    flags = {'CFLAGS', '-std=gnu99 -fopenmp -O3'; 'LDFLAGS', '-fopenmp'};
end

mexArgs = {'-largeArrayDims', ['-I' srcDir], ['-I' generatedDir]};

% For Octave the flags are set in the environment, while MATLAB parses them
% as arguments
for flag = 1:size(flags, 1)
    if strcmp(env, 'OCTAVE')
        preFlagContent = evalc(['mkoctfile -p ' flags{flag, 1}]);

        if ~isempty(preFlagContent)
            preFlagContent = preFlagContent(1:end - 1); % Strip newline
        end

        newContent = [preFlagContent ' ' flags{flag, 2}];
        setenv(flags{flag, 1}, newContent);
        matRad_cfg.dispDebug('Set compiler flag %s to %s\n', flags{flag, 1}, newContent);
    else
        mexArgs{end + 1} = sprintf('%s="$%s %s"', flags{flag, 1}, flags{flag, 1}, flags{flag, 2}); %#ok<AGROW>
    end
end

mexArgs = [mexArgs, sourceFiles];

outDir = fullfile(tempname, 'ompMC-mex');
mkdir(outDir);

matRad_cfg.dispInfo('Building ompMC with mex...\n');
matRad_cfg.dispDebug('Compiler call: mex %s\n', strjoin(mexArgs, ' '));

% mex writes to the working directory, and -outdir is not understood by both
% environments
currDir = pwd;
returnToCurrDir = onCleanup(@() cd(currDir));
cd(outDir);

mex(mexArgs{:});

clear returnToCurrDir;

mexFile = matRad_ompMCfindMexFile(outDir);

if isempty(mexFile)
    matRad_cfg.dispError('The mex call produced no omc_matrad mex file in %s.', outDir);
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCwriteVersionHeader(sourceFolder, generatedDir)
% Fills in src/omc_version.h.in the way CMake would, taking the version from
% the project() call that ompMC keeps it in

matRad_cfg = MatRad_Config.instance();

cmakeLists = fileread(fullfile(sourceFolder, 'CMakeLists.txt'));
versionToken = regexp(cmakeLists, 'project\s*\(\s*ompMC\s+VERSION\s+(\d+)\.(\d+)\.(\d+)', 'tokens', 'once');

if isempty(versionToken)
    matRad_cfg.dispError('Could not read the ompMC version from %s.', fullfile(sourceFolder, 'CMakeLists.txt'));
end

template = fileread(fullfile(sourceFolder, 'src', 'omc_version.h.in'));

replacements = { ...
                '@PROJECT_VERSION_MAJOR@', versionToken{1}; ...
                '@PROJECT_VERSION_MINOR@', versionToken{2}; ...
                '@PROJECT_VERSION_PATCH@', versionToken{3}; ...
                '@PROJECT_VERSION@',       strjoin(versionToken, '.')};

for i = 1:size(replacements, 1)
    template = strrep(template, replacements{i, 1}, replacements{i, 2});
end

fid = fopen(fullfile(generatedDir, 'omc_version.h'), 'w');

if fid < 0
    matRad_cfg.dispError('Could not write %s.', fullfile(generatedDir, 'omc_version.h'));
end

fprintf(fid, '%s', template);
fclose(fid);
end

%% ------------------------------------------------------------------------
function mexFile = matRad_ompMCfindMexFile(folder)
% The extension differs per platform and environment, so the file is looked
% up by name rather than constructed

mexFile = '';

candidates = dir(fullfile(folder, ['omc_matrad.' mexext]));

if ~isempty(candidates)
    mexFile = fullfile(folder, candidates(1).name);
end
end

%% ------------------------------------------------------------------------
function tf = matRad_ompMCisLoaded()

try
    tf = any(strcmp(inmem('-mex'), 'omc_matrad'));
catch
    tf = false;
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCmoveContent(fromDir, toDir)
% Moves everything in fromDir into toDir, which is how an unpacked release
% becomes an installation. Both are on the same volume by construction: on
% Windows octave moves a directory by asking the system to, which it refuses
% to do across drives.

matRad_cfg = MatRad_Config.instance();

if ~exist(toDir, 'dir')
    mkdir(toDir);
end

entries = dir(fromDir);

for i = 1:numel(entries)
    if any(strcmp(entries(i).name, {'.', '..'}))
        continue
    end

    [success, msg] = movefile(fullfile(fromDir, entries(i).name), fullfile(toDir, entries(i).name));

    if ~success
        matRad_cfg.dispError('Could not move %s into %s: %s', entries(i).name, toDir, msg);
    end
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCclean(installDir)
% Removes an installation, keeping the installer itself - it is the only
% tracked file in this folder

if ~exist(installDir, 'dir')
    return
end

entries = dir(installDir);

for i = 1:numel(entries)
    if any(strcmp(entries(i).name, {'.', '..'})) || ~isempty(regexp(entries(i).name, '\.m$', 'once'))
        continue
    end

    entryPath = fullfile(installDir, entries(i).name);

    if entries(i).isdir
        % matRad adds everything under thirdParty to the path on startup, and
        % octave complains about every folder that disappears underneath it
        warnState = warning();
        warning('off', 'all');
        rmpath(entryPath);
        warning(warnState);

        matRad_ompMCrmdir(entryPath);
    else
        delete(entryPath);
    end
end

% Windows does not let go of a mex file another session has loaded, and
% installing over it would fail later with whatever the copy makes of that
leftoverMex = matRad_ompMCfindMexFile(fullfile(installDir, 'bin'));

if ~isempty(leftoverMex)
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError(['%s could not be removed. It is most likely loaded by another ' ...
                          'MATLAB or Octave session - close those and try again.'], leftoverMex);
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCrmdir(folder)

if exist(folder, 'dir')
    [status, msg] = rmdir(folder, 's');

    if status ~= 1
        matRad_cfg = MatRad_Config.instance();
        matRad_cfg.dispWarning('Could not remove %s: %s', folder, msg);
    end
end
end

%% ------------------------------------------------------------------------
function matRad_ompMCwriteOrigin(installDir, release, how)
% Records where this installation came from, next to the license it came
% under - the folder is not tracked, so nothing else says so

fid = fopen(fullfile(installDir, 'ORIGIN.txt'), 'w');

if fid < 0
    return
end

fprintf(fid, 'ompMC %s\n', release.tag);
fprintf(fid, '%s\n\n', release.repoUrl);
fprintf(fid, 'Installed by matRad_installOmpMC as %s.\n\n', how);
fprintf(fid, 'ompMC is licensed under the GNU General Public License v3, see the LICENSE\n');
fprintf(fid, 'file next to this one. It is not part of matRad, which is licensed under the\n');
fprintf(fid, '3-clause BSD license, and it is not distributed with matRad - this folder was\n');
fprintf(fid, 'populated on this machine and is not tracked by git.\n');
fclose(fid);
end
