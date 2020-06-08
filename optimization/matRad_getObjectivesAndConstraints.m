function classNames = matRad_getObjectivesAndConstraints()
% matRad steering information generation
% 
% call
%   classNames = matRad_getObjectivesAndConstraints()
%
%
% output
%   classNames:     contains class names (row 1) and display descriptions 
%                   (row 2) of all available objectives
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

env = matRad_getEnvironment();

if strcmp(env,'MATLAB') %use the package methodology for Matlab (stable)
    mpkgObjectives = meta.package.fromName('DoseObjectives');
    mpkgConstraints = meta.package.fromName('DoseConstraints');
    classList = [mpkgObjectives.ClassList; mpkgConstraints.ClassList];

    classList = classList(not([classList.Abstract]));

    %Now get the "internal" name from the properties
    classNames = cell(2,numel(classList));
    for clIx = 1:numel(classList)
        cl = classList(clIx);
        pList = cl.PropertyList; %Get List of all properties
        pNameIx = arrayfun(@(p) strcmp(p.Name,'name'),pList); %get index of the "name" property
        p = pList(pNameIx); %select name property
        pName = p.DefaultValue; % get value / name
        classNames(:,clIx) = {cl.Name; pName}; %Store class name and display name
    end
else
    currDir = fileparts(mfilename('fullpath'));
    objectiveDir = [currDir filesep '+DoseObjectives'];
    constraintDir = [currDir filesep '+DoseConstraints'];
       
    objFiles = dir(objectiveDir);
    for i=1:numel(objFiles)
      objFiles(i).pkgName = 'DoseObjectives';
    end
    constrFiles = dir(constraintDir);
    for i=1:numel(objFiles)
      constrFiles(i).pkgName = 'DoseConstraints';
    end
    
    allFiles = [objFiles; constrFiles];
   
    [defNames,dispNames] = arrayfun(@testForClass,allFiles,'UniformOutput',false);
    
    classNames(1,:) = defNames;
    classNames(2,:) = dispNames;
    
    selectIx = cellfun(@isempty,classNames);
    classNames = classNames(:,~selectIx(1,:));    
end
end

function [defName,dispName] = testForClass(potentialClassFile)
    
    
    try
%        fPath = potentialClassFile.folder;
        [~,fName,~] = fileparts(potentialClassFile.name);
        fType = potentialClassFile.pkgName;
      
        clTmp = meta.class.fromName([fType '.' fName]);
        defName = clTmp.Name;
       
        pList = clTmp.PropertyList;
        
        pNameIx = cellfun(@(p) strcmp(p.Name,'name'),pList);
        p = pList(pNameIx);
        dispName = p{1}.DefaultValue                ;
        
        
    catch
        defName = '';
        dispName = '';
    end      
end

