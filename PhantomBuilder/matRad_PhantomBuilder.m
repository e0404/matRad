classdef matRad_PhantomBuilder < handle
    % matRad_PhantomBuilder
    % Class that helps to create radiotherapy phantoms 
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2023 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        volumes = {};
        ct;
        cst = {};    
    end

    methods
        function obj = matRad_PhantomBuilder(ct)
            obj.ct = ct;
        end

        %functions to create Targets %TODO: Option to extend volumes
        function addCubicTarget(obj,name,dimensions,varargin)
            obj.volumes(end+1) = {matRad_CubicVOI(name,'TARGET',dimensions,varargin{:})}
        end

        function addSphericalTarget(obj,name,radius,varargin)
            obj.volumes(end+1) = {matRad_SphericalVOI(name,'TARGET',radius,varargin{:})} 
        end


        
        function addCubicOAR(obj,name,dimensions,varargin)
            obj.volumes(end+1) = {matRad_CubicVOI(name,'OAR',dimensions,varargin{:})}
        end

        function addSphericalOAR(obj,name,radius,varargin)
            obj.volumes(end+1) ={matRad_SphericalVOI(name,'OAR',radius,varargin{:})}
        end
        

        function updatecst(obj) %after initializing the volumes used to create the cst struct
            for i = 1:numel(obj.volumes)
                obj.cst = obj.volumes{i}.initializeParameters(obj.ct,obj.cst)    
            end
        end

        function updatect(obj) %update the ct data %TODO: What if different HU unites than water are required?
            for i = 1:size(obj.cst,1)
                vIxVOI = obj.cst{i,4}{1};
                obj.ct.cubeHU{1}(vIxVOI) = 0; % assign HU of water
            end
        end

    end    
end