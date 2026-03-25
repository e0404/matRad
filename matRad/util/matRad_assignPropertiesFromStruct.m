function matRad_assignPropertiesFromStruct(obj,propertiesStruct,overwrite,fieldChangedWarningMessage)
% matRad helper function to configure object from structure
%
% call:
%   obj = matRad_assignPropertiesFromStruct(obj,propertiesStruct,overwrite,fieldChangedWarningMessage)
%   obj = matRad_assignPropertiesFromStruct(obj,propertiesStruct,overwrite)
%   obj = matRad_assignPropertiesFromStruct(obj,propertiesStruct)
%
% input:
%   obj:                        Object to be configured
%   propertiesStruct:           Structure containing properties to be
%                               assigned
%   overwrite:                  (optional) Boolean flag to overwrite
%                               existing properties (default: true)
%   fieldChangedWarningMessage: (optional) Custom warning message for
%                               changed fields
%
% output:
%   resultGUI:  struct containing optimized fluence vector, dose, and (for
%               biological optimization) RBE-weighted dose etc.
%   optimizer:  Used Optimizer Object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team.
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

if nargin < 4
    fieldChangedWarningMessage = '';
end

if nargin < 3
    overwrite = true;
end

if matRad_cfg.isOctave
    c2sWarningState = warning('off','Octave:classdef-to-struct');
end

fields = fieldnames(propertiesStruct);

for i = 1:numel(fields)
    field = fields{i};
    try
        if matRad_ispropCompat(obj,field)
            obj.(field) = matRad_recursiveFieldAssignment(obj.(field),propertiesStruct.(field),overwrite,fieldChangedWarningMessage,field);
        else
            matRad_cfg.dispWarning('Tried to assign nonexisting property ''%s'' from property struct to Object of class %s!',field,class(obj));
        end
    catch ME
        % catch exceptions when the class has no properties defined in the
        % struct.
        switch ME.identifier
            case {'MATLAB:class:noPublicField','MATLAB:class:noSetMethod'}
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Could not assign field %s to object of class %s: %s', field, class(obj), ME.message);
            otherwise
                matRad_cfg.dispWarning('Problem while setting up object of class %s from struct:%s %s',class(obj),field,ME.message);
        end
    end
end

if matRad_cfg.isOctave
    warning(c2sWarningState.state,'Octave:classdef-to-struct');
end
