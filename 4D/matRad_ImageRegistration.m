classdef (Abstract) matRad_ImageRegistration
    % matRad_imageRegistration. Superclass for objectives and constraints
    % This is the superclass for all objectives and constraints to enable easy
    % one-line identification.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties (Abstract, Constant)
        name                %Display name of the Objective. Needs to be implemented in sub-classes.
    end
    
    properties (Abstract)
        refScen
        metadata
    end
    
    methods
        function obj = matRad_ImageRegistration(dataStruct)
            if nargin > 0 && ~isempty(dataStruct) && isstruct(dataStruct)
                obj = assignCommonPropertiesFromStruct(obj,dataStruct);
            end
        end
        
        % Overload the struct function to return a struct with general
        % the objective / constraint
        function s = struct(obj)
            s.className = class(obj);
            s.ct = obj.ct;
            s.cst = obj.cst;
            s.refScen = obj.refScen;
            s.metadata = obj.metadata;
        end
    end
    
	methods (Abstract)
        calcDVF(obj,ct,cst)
        propContours(ct,cst)
	end
       
end

