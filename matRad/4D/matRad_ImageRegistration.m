classdef (Abstract) matRad_ImageRegistration
    % matRad_ImageRegistration Abstract superclass for 4D image registration
    %   Defines the common interface for image registration classes that compute
    %   deformation vector fields and propagate contours between CT scenarios.
    %
    % References
    %   -
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

    properties (Abstract, Constant)
        name
    end

    properties (Abstract)
        ct
        cst
        refScen
        metadata
    end

    methods

        function obj = matRad_ImageRegistration(dataStruct)
            if nargin > 0 && ~isempty(dataStruct) && isstruct(dataStruct)
                obj = obj.assignCommonPropertiesFromStruct(dataStruct);
            end
        end

        function s = struct(obj)
            s = struct();
            s.className = class(obj);
            s.ct = obj.ct;
            s.cst = obj.cst;
            s.refScen = obj.refScen;
            s.metadata = obj.metadata;
        end

    end

    methods (Access = private)

        function obj = assignCommonPropertiesFromStruct(obj, s)
            names = fieldnames(s);
            for i = 1:numel(names)
                try
                    obj.(names{i}) = s.(names{i});
                catch
                end
            end
        end

    end

    methods (Abstract)
        ct = calcDVF(obj)
        [ct, cst] = propContours(obj)
    end
end
