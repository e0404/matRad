function [classList] = matRad_identifyClassesByConstantProperties(metaClasses,primaryPropertyName,varargin)
% matRad_identifyClassesByProperty: Helper function to identify classes by
% property
%   This method identifies classes based on a primary property and optional
%   additional properties.
%
% call:
%   classList = matRad_identifyClassesByProperty(metaClasses,
%   primaryPropertyName) classList =
%       matRad_identifyClassesByProperty(metaClasses, primaryPropertyName,
%       'defaults', {defaultClasses})
%   classList =
%       matRad_identifyClassesByProperty(metaClasses, primaryPropertyName,
%       'additionalPropertyNames', {additionalProperties})
%
% inputs:
%   - metaClasses: A cell array of meta.class objects representing the
%     classes to be identified. 
%   - primaryPropertyName: The name of the primary property used for 
%     identification.
%
% optional Parameter Inputs:
%   - defaults: A cell array of default classes that should be listed
%     first. 
%   - additionalPropertyNames: A cell array of additional property names
%     used for identification.
%
% outputs:
%   - classList: A structure array containing the identified classes.
%       - primaryPropertyName: The values of the primary property for each
%       class. - additionalPropertyNames: The values of the additional
%       properties for each class. - className: The names of the identified
%       classes. - handle: The constructor handles of the identified
%       classes.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance(); % Get the instance of the MatRad_Config class

p = inputParser; % Create an input parser object

% Define the required inputs for the input parser
p.addRequired('metaClasses',@(x) iscell(x) && all(cellfun(@(x) isa(x,'meta.class'),x)));
p.addRequired('primaryPropertyName',@(x) ischar(x) && isrow(x) && ~isempty(x));

% Define the optional inputs for the input parser
p.addParameter('defaults',{},@(x) iscell(x) && all(cellfun(@(x) ischar(x) && isrow(x),x)));
p.addParameter('additionalPropertyNames',{},@(x) iscell(x) && all(cellfun(@(x) ischar(x) && isrow(x),x)));

% Parse the input arguments
p.parse(metaClasses,primaryPropertyName,varargin{:});
defaults = p.Results.defaults; 
additionalPropertyNames = p.Results.additionalPropertyNames; 

% Get the values of the primary property for each class
primaryPropertyValueList = getPropertyValues(metaClasses,primaryPropertyName,matRad_cfg.isMatlab); 

% Get the values of the additional properties for each class
additionalPropertyValueList = cell(numel(additionalPropertyNames),numel(metaClasses));
for i = 1:numel(additionalPropertyNames)
    additionalPropertyName = additionalPropertyNames{i};
    additionalPropertyValueList(i,:) = getPropertyValues(metaClasses,additionalPropertyName,matRad_cfg.isMatlab);
end

% Get the names of the identified classes and their constructor handles
classNameList = cellfun(@(mc) mc.Name,metaClasses,'UniformOutput',false); 
constructorHandleList = cellfun(@(namestr) str2func(namestr),classNameList,'UniformOutput',false); 

% Create a default sort pattern
sortPattern = 1:numel(primaryPropertyValueList); 

% Make sure the default engines are the first ones listed
if ~isempty(defaults)
    findDefaultIx = [];
    for defaultClassIx = 1:length(defaults)
        defaultClass = defaults{defaultClassIx};
        foundIx = find(strcmp(defaultClass,primaryPropertyValueList));
        findDefaultIx(end+1:end+numel(foundIx)) = foundIx;
    end

    if ~isempty(findDefaultIx)
        sortPattern = [sortPattern(findDefaultIx), sortPattern];
        sortPattern = unique(sortPattern,'stable');
    end
end

% Sort the property values, class names, and constructor handles
primaryPropertyValueList = primaryPropertyValueList(sortPattern); 
additionalPropertyValueList = additionalPropertyValueList(:,sortPattern);
classNameList = classNameList(sortPattern); 
constructorHandleList = constructorHandleList(sortPattern);

% Create a structure array containing the identified classes
classList = cell2struct([primaryPropertyValueList; additionalPropertyValueList; classNameList; constructorHandleList],{primaryPropertyName,additionalPropertyNames{:},'className','handle'});

end

%Helper function to get the property values dependent on the used environment
function valueList = getPropertyValues(metaClasses,propertyName,isMatlab)
    if isMatlab
        % Get the default values of the property for each class in MATLAB
        valueList = cellfun(@(mc) mc.PropertyList(strcmp({mc.PropertyList.Name}, propertyName)).DefaultValue,metaClasses,'UniformOutput',false); 
    else
        % Get the default values of the property for each class in Octave
        valueList = cellfun(@(mc) mc.PropertyList{find(cellfun(@(p) strcmp(p.Name, propertyName),mc.PropertyList))}.DefaultValue,metaClasses,'UniformOutput',false);
    end
end
