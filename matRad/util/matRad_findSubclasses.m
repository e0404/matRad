function classList = matRad_findSubclasses(superClass,varargin)
% matRad_findSubclasses: Helper function to find subclasses
%   This method can find subclasses to a class within package folders and
%   normal folders. It is not very fast, but Matlab & Octave compatible, so
%   it is advised to cache classes once scanned.
%
% call:
%   classList = matRad_findSubclasses(superClass);
%   classList =
%       matRad_findSubclasses(superClass,'package',{packageNames})
%   classList =
%       matRad_findSubclasses(superClass,'folders',{folderNames})
%   classList =
%       matRad_findSubclasses(superClass,'includeAbstract',true/false)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addRequired('superClass',@(x) ischar(x) || isa(x,'meta.class'));
p.addParameter('packages',{},@(x) iscell(x) && (iscellstr(x) || all(cellfun(@(y) isa(y,'meta.package'),x))));
p.addParameter('folders',{},@(x) iscellstr(x) && all(isfolder(x)));
p.addParameter('includeAbstract',false,@(x) isscalar(x) && islogical(x));
p.addParameter('includeSubfolders',false,@(x) isscalar(x) && islogical(x));
p.addParameter('usePath',false,@(x) islogical(x) && isscalar(x));

p.parse(superClass,varargin{:});

superClass = p.Results.superClass;
packages = p.Results.packages;
folders = p.Results.folders;
includeAbstract = p.Results.includeAbstract;
includeSubfolders = p.Results.includeSubfolders;
usePath = p.Results.usePath;

%Create MetaObject if char was given
if ischar(superClass)
    superClass = meta.class.fromName(superClass);
end

%Check if path is used
if usePath
    addPath = [path pathsep];
else
    addPath = '';
end

for folderIx = 1:length(folders)
    %Generate subfolder path or only add current folder to class search
    %path
    if includeSubfolders
        addPath = [addPath genpath(folders{folderIx}) pathsep];
    else
        addPath = [addPath folders{folderIx} pathsep];
    end
end
folders = strsplit(addPath,pathsep);

if all(cellfun(@(y) isa(y,'meta.package'),packages))
    packages = cellfun(@(y) y.Name,packages);
end

matRad_cfg = MatRad_Config.instance();

classList = {};

%%Collect package classes
for packageIx = 1:length(packages)
    package = packages{packageIx};
    if matRad_cfg.isMatlab
        mp = meta.package.fromName(package);
        classList = [classList; num2cell(mp.ClassList)];
    elseif matRad_cfg.isOctave
        p = strsplit(path,pathsep);
        ix = find(cellfun(@(f) isfolder([f filesep '+' package]),p));
        if ~isscalar(ix)
            matRad_cfg.dispError('Could not uniquely identify package folder ''+%s'' (%d folders found with that name!)',package,numel(ix));
        end
        classList = [classList; matRad_getClassesFromFolder(p{ix},package)];
    else
        matRad_cfg.dispError('Environment %s unknown!',matRad_cfg.env);
    end
end

%Collect classes from folders
for folderIx = 1:length(folders)
    folder = folders{folderIx};
    classList = [classList; matRad_getClassesFromFolder(folder)];
end

%Throw out abstract classes if requested
if ~includeAbstract
    isNotAbstract = cellfun(@(mc) ~mc.Abstract,classList);
    classList = classList(isNotAbstract);
end

%Now test if it is a subclass

inherits = cellfun(@(mc) matRad_checkInheritance(mc,superClass),classList);
classList = classList(inherits);

%We want classes to be a row vector, but the meta.class lists column
%vector
if ~isrow(classList)
    classList = classList';
end
end

function metaClassList = matRad_getClassesFromFolder(folder,packageName)
% matRad_getClassesFromFolder: collects classes from a folder. packageName
% can be passed for Octave compatibility
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Sanity check to avoid octave calling what when argument empty
if isempty(folder)
    metaClassList = {};
    return;
end

%We accept the specific package name for Octave compatability
if nargin < 2
    packageNameWithDot= '';
    fullFolder = folder;
else
    packageNameWithDot = [packageName '.'];
    fullFolder =  [folder filesep '+' packageName];
end

folderInfo = what(fullFolder);
if ~isempty(folderInfo)
    [~,potentialClasses] = cellfun(@fileparts,[folderInfo.m(:)',folderInfo.p(:)'],'UniformOutput',false); %Potential class files
    
    %In octave the what function returns the class folders with an '@'
    classFolders = folderInfo.classes;
    classFolders = cellfun(@(f) erase(f,'@'),classFolders,'UniformOutput',false);

    potentialClasses = [potentialClasses, classFolders];

else
    potentialClasses = {};
end

metaClassList = {};

for potentialClassIx = 1:length(potentialClasses)
    potentialClass = potentialClasses{potentialClassIx};
    % get meta class with package name infront
    mc = meta.class.fromName([packageNameWithDot potentialClass]);

    % Check if we found a class
    if isempty(mc)
        continue;
    end

    metaClassList{end+1} = mc;
end
metaClassList = transpose(metaClassList);
end

function inherits = matRad_checkInheritance(metaClass,superClass)
% matRad_checkInheritance: checks class inheritance of a superclass
% recursively. Would be easier with superclass method, but this is not
% available in Octave
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

if matRad_cfg.isMatlab
    scs = superclasses(metaClass.Name);
    inherits = any(strcmp(scs,superClass.Name));
else
    superClasses = metaClass.SuperclassList;

    if isempty(superClasses)
        inherits = false;
    else
        inherits = any(cellfun(@(x) strcmp(x.Name,superClass.Name),superClasses));
        if inherits == false
            for ix = 1:numel(superClasses)
                inherits = matRad_checkInheritance(superClasses{ix},superClass);
            end
        end
    end
end

end
