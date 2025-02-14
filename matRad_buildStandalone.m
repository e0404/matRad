function buildResult = matRad_buildStandalone(varargin)
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

%Setup Input Parsing
p = inputParser;
p.addParameter('isRelease',false,@islogical);                                           %By default we compile a snapshot of the current branch
p.addParameter('compileWithRT',false,@islogical);                                       %By default we don't package installers with runtime
p.addParameter('buildDir',fullfile(matRadRoot,'build'),@(x) ischar(x) || isstring(x));  %Build directory
p.addParameter('verbose',false,@islogical);
p.addParameter('docker',false,@islogical);
p.addParameter('java',false,@islogical);
p.addParameter('python',false,@islogical);
p.addParameter('json',[],@(x) isempty(x) || ischar(x) || isstring(x));

%Parse and manage inputs
p.parse(varargin{:});
isRelease = p.Results.isRelease;
compileWithRT = p.Results.compileWithRT;
buildDir = p.Results.buildDir;
buildDocker = all(p.Results.docker);
buildJava   = all(p.Results.java);
buildPython = all(p.Results.python);
json  = p.Results.json;

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

%Display OS for debugging and information
arch = computer('arch');
fprintf('Build Architecture: %s\n',arch);
archcheck = string([ispc,isunix,ismac]).cellstr();
fprintf('pc\t\tunix\tmac\n%s\t%s\t%s\n',archcheck{:});

%Check build directory
try
    mkdir(buildDir);
catch ME
    error(ME.identifier,'Could not create build directory %s\n Error:',buildDir,ME.message);
end

%Docker build?
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

buildResult.version = vernumInstall;
buildResult.versionDetail = versionFull;

%% Set Options and Compile
try
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
        resultsStandalone = compiler.build.standaloneWindowsApplication(buildOpts);
    else
        resultsStandalone = compiler.build.standaloneApplication(buildOpts);
    end

    buildResult.standalone.compiledFiles = resultsStandalone.Files;

catch ME
    warning(ME.identifier,'Failed to compile standalone due to %s',ME.message);

end

%% Package
if ispc
    readmeFile = 'readme_standalone_windows.txt';
elseif ismac
    readmeFile = 'readme_standalone_mac.txt';
else
    readmeFile = 'readme_standalone_linux.txt';
end
installerId = arch;

try
    packageOpts = compiler.package.InstallerOptions(resultsStandalone,...
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
        'InstallerName',sprintf('matRad_installer_%s_v%s',installerId,vernumApp),...
        'OutputDir',fullfile(buildDir,'installer'),...
        'RuntimeDelivery',rtOption,...
        'Summary','matRad is an open source treatment planning system for radiation therapy written in Matlab.',...
        'Version',vernumInstall,...
        'Verbose',verbose);
    compiler.package.installer(resultsStandalone,'Options',packageOpts);

    outFiles = dir([packageOpts.OutputDir filesep packageOpts.InstallerName '.*']);
    outFiles = arrayfun(@(f) fullfile(f.folder,f.name),outFiles,'UniformOutput',false);

    buildResult.standalone.installerFiles = outFiles;
catch ME
    warning(ME.identifier,'Failed to package standalone installer due to %s!',ME.message);
end

%% Python
if buildPython
    functionFiles = {'matRadGUI.m',...
        'matRad/matRad_generateStf.m',...
        'matRad/matRad_calcDoseForward.m',...
        'matRad/matRad_calcDoseInfluence.m',...
        'matRad/matRad_fluenceOptimization.m',...
        'matRad/matRad_directApertureOptimization.m',...
        'matRad/matRad_sequencing.m',...
        'matRad/matRad_planAnalysis.m'};

    sampleGenerationFiles = {'matRad.m'};

    try
        pythonOpts = compiler.build.PythonPackageOptions(functionFiles,...
            'AdditionalFiles',{'matRad','thirdParty','matRad_rc.m'},...
            'SampleGenerationFiles',sampleGenerationFiles,...
            'Verbose',verbose, ...
            'PackageName','pyMatRad',...
            'OutputDir',fullfile(buildDir,'python'));

        resultsPython = compiler.build.pythonPackage(pythonOpts);

        buildResult.python.packageFiles = resultsPython.Files;
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
        resultsJava = compiler.build.javaPackage(javaOpts);
        buildResult.java.packageFiles = resultsJava.Files;
    catch ME
        warning(ME.identifier,'Java build failed due to %s!',ME.message);
    end
end

%% Docker
if buildDocker
    try
        if isRelease
            imageName = ['matRad:' vernumInstall];
        else
            imageName = 'matRad:develop';
        end
        dockerOpts = compiler.package.DockerOptions(results,...
            'AdditionalInstructions','',...
            'AdditionalPackages','',...
            'ContainerUser','appuser',...
            'DockerContext',fullfile(buildDir,'docker'),...
            'ExecuteDockerBuild','on',...
            'ImageName',imageName);

        compiler.package.docker(results,'Options',dockerOpts);
        buildResult.docker.image = dockerOpts.ImageName;
    catch ME
        warning(ME.identifier,'Java build failed due to %s!',ME.message);
    end
end

if ~isempty(json)
    try
        fH = fopen(json,'w');
        jsonStr = jsonencode(buildResult,"PrettyPrint",true);
        fwrite(fH,jsonStr);
        fclose(fH);
    catch ME
        warning(ME.identifier,'Could not open JSON file for writing: %s',ME.message);
    end
end
