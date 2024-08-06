function matRad_buildStandalone(varargin)
% Compiles the standalone exectuable & packages installer using Matlab's
% Compiler Toolbox
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = matRad_rc;

matRadRoot = matRad_cfg.matRadRoot;
standaloneFolder = fullfile(matRadRoot,'standalone');

p = inputParser;

p.addParameter('isRelease',false,@islogical);                                           %By default we compile a snapshot of the current branch
p.addParameter('compileWithRT',false,@islogical);                                       %By default we don't package installers with runtime
p.addParameter('buildDir',fullfile(matRadRoot,'build'),@(x) ischar(x) || isstring(x));  %Build directory
p.addParameter('verbose',false,@islogical);
p.addParameter('docker',false,@islogical);
p.addParameter('java',false,@islogical);
p.addParameter('python',false,@islogical);

p.parse(varargin{:});

isRelease = p.Results.isRelease;
compileWithRT = p.Results.compileWithRT;
buildDir = p.Results.buildDir;
buildDocker = all(p.Results.docker);
buildJava   = all(p.Results.java);
buildPython = all(p.Results.python);

if all(p.Results.verbose)
    verbose = 'On';
else
    verbose = 'Off';
end

if compileWithRT
    rtOption = 'installer';
else
    rtOption = 'web';
end

try 
    mkdir(buildDir);
catch ME
    error(ME.identifier,'Could not create build directory %s\n Error:',buildDir,ME.message);
end

if buildDocker && isunix && ~ismac
    warning('Can''t build docker container. Only works on linux!');
    buildDocker = false;
end


[~,versionFull] = matRad_version();

if isRelease
    vernumApp = sprintf('%d.%d.%d',versionFull.major,versionFull.minor,versionFull.patch);
    vernumInstall = vernumApp;
    %vernumInstall = sprintf('%d.%d',versionFull.major,versionFull.minor);
else
    vernumApp = sprintf('%d.%d.%d.65534',versionFull.major,versionFull.minor,versionFull.patch);
    vernumInstall = vernumApp;
    %vernumInstall = sprintf('%d.%d',versionFull.major,versionFull.minor);
end


%elseif isempty(versionFull.commitID)
%    vernum = sprintf('%d.%d.%d.dev',versionFull.major,versionFull.minor,versionFull.patch);
%else
%    vernum = sprintf('%d.%d.%d.%s/%s',versionFull.major,versionFull.minor,versionFull.patch,versionFull.branch,versionFull.commitID(1:8));
%end

%% Set Options
buildOpts = compiler.build.StandaloneApplicationOptions('matRadGUI.m',...
    'OutputDir',buildDir,...
    'ExecutableIcon',fullfile(standaloneFolder,'matRad_icon.png'),...
    'AutodetectDataFiles','on',...
    'AdditionalFiles',{'matRad','thirdParty','matRad_rc.m'},...
    'EmbedArchive','on',...
    'ExecutableName','matRad',...
    'ExecutableSplashScreen',fullfile(standaloneFolder,'matRad_splashscreen.png'),...
    'ExecutableVersion',vernumApp,...
    'TreatInputsAsNumeric','off',...
    'Verbose',verbose);

if ispc
    results = compiler.build.standaloneWindowsApplication(buildOpts);
else
    results = compiler.build.standaloneApplication(buildOpts);
end

%% Package
if ispc
    readmeFile = 'readme_standalone_windows.txt';
    installerId = 'Win64';
elseif ismac
    readmeFile = 'readme_standalone_mac.txt';
    [~,result] = system('uname -m');
    if any(strfind(result,'ARM64')) %is m1mac
        installerId = 'MacARM';
    else
        installerId = 'Mac64';
    end

else
    readmeFile = 'readme_standalone_linux.txt';
    installerId = 'Linux64';
end

packageOpts = compiler.package.InstallerOptions(results,...
    'AdditionalFiles',{fullfile(matRadRoot,'AUTHORS.txt'),fullfile(matRadRoot,'LICENSE.md'),fullfile(standaloneFolder,readmeFile)},...
    'ApplicationName','matRad',...
    'AuthorCompany','German Cancer Research Center (DKFZ)',...
    'AuthorEmail','contact@matRad.org',...
    'AuthorName','matRad development team @ DKFZ',...
    'Description','matRad is an open source software for radiation treatment planning of intensity-modulated photon, proton, and carbon ion therapy started in 2015 in the research group "Radiotherapy Optimization" within the Department of Medical Physics in Radiation Oncology at the German Cancer Research Center - DKFZ.', ... %\n\nmatRad targets education and research in radiotherapy treatment planning, where the software landscape is dominated by proprietary medical software. As of August 2022, matRad had more than 130 forks on GitHub and its development paper was cited more than 160 times (according to Google Scholar). matRad is entirely written in MATLAB and mostly compatible to GNU Octave.',...
    'InstallationNotes','matRad contains precompiled libraries and thirdParty software. In some cases, those precompiled libraries may not run out of the box on your system. Please contact us should this be the case. For the Third-Party licenses, check the thirdParty subfolder in the matRad installation directory.',...
    'InstallerIcon',fullfile(standaloneFolder,'matRad_icon.png'),...
    'InstallerLogo',fullfile(standaloneFolder,'matRad_installscreen.png'),...
    'InstallerSplash',fullfile(standaloneFolder,'matRad_splashscreen.png'),...
    'InstallerName',sprintf('matRad_installer%s_v%s',installerId,vernumApp),...
    'OutputDir',fullfile(buildDir,'installer'),...
    'RuntimeDelivery',rtOption,...
    'Summary','matRad is an open source treatment planning system for radiation therapy written in Matlab.',...
    'Version',vernumInstall,...
    'Verbose',verbose);

compiler.package.installer(results,'Options',packageOpts);

%% Python
if buildPython
    pythonOpts = compiler.build.PythonPackageOptions('matRadGUI.m',...
        'AdditionalFiles',{'matRad','thirdParty','matRad_rc.m'},...
        'Verbose',verbose, ...
        'PackageName','pyMatRad',...
        'OutputDir',fullfile(buildDir,'python'));
    try
        results = compiler.build.pythonPackage(pythonOpts);
    catch ME
        warning(ME.identifier,'Java build failed due to %s!',ME.message);
    end
end

%% Java
if buildJava
    javaOpts = compiler.build.JavaPackageOptions('matRadGUI.m',...
        'AdditionalFiles',{'matRad','thirdParty','matRad_rc.m'},...
        'Verbose',verbose);
    try
        results = compiler.build.javaPackage(javaOpts);
    catch ME
        warning(ME.identifier,'Java build failed due to %s!',ME.message);
    end
end

%% Docker
if ~buildDocker
    return;
end

dockerOpts = compiler.package.DockerOptions(results,...
    'AdditionalInstructions','',...
    'AdditionalPackages','',...
    'ContainerUser','appuser',...
    'DockerContext',fullfile(buildDir,'docker'),...
    'ExecuteDockerBuild','on',...
    'ImageName','matRad',...
    'RuntimeImage','mcrimage');

compiler.package.docker(results,'Options',dockerOpts);


   








