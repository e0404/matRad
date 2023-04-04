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

    properties (Access = public)
        volumes = {};
    end

    properties (Access = private) 
        ct;
        cst = {};    
    end

    methods (Access = public)
        function obj = matRad_PhantomBuilder(ctDim,ctResolution,numOfCtScen)
            obj.ct = {};
            obj.ct.cubeDim = ctDim;
            obj.ct.resolution.x = ctResolution(2);
            obj.ct.resolution.y = ctResolution(1);
            obj.ct.resolution.z = ctResolution(3);
            obj.ct.numOfCtScen = numOfCtScen;
            obj.ct.cubeHU{1} = ones(obj.ct.cubeDim) * -1000;
        end

        %functions to create Targets %TODO: Option to extend volumes
        function addBoxTarget(obj,name,dimensions,varargin)
            % Adds a box target
            %
            % input:
            %   name:           Name of VOI as string
            %   dimensions:     Dimensions of the box as [x,y,z] array
            %
            % Name-Value pairs:
            %   'offset':       The offset of the VOI with respect to the
            %                   center of the geometry as [x,y,z] array
            %   'objectives':   Either a single objective or a cell array
            %                   of objectives
            %   'HU':           Houndsfield unit of the volume

            obj.volumes(end+1) = {matRad_PhantomVOIBox(name,'TARGET',dimensions,varargin{:})};
            obj.updatecst();
        end

        function addSphericalTarget(obj,name,radius,varargin)
            % Adds a spherical target
            %
            % input:
            %   name:           Name of VOI as string
            %   radius:         Radius of the sphere
            %
            % Name-Value pairs:
            %   'offset':       The offset of the VOI with respect to the
            %                   center of the geometry as [x,y,z] array
            %   'objectives':   Either a single objective or a cell array
            %                   of objectives
            %   'HU':           Houndsfield unit of the volume

            obj.volumes(end+1) = {matRad_PhantomVOISphere(name,'TARGET',radius,varargin{:})};
            obj.updatecst();
        end

        
        function addBoxOAR(obj,name,dimensions,varargin)
            % Adds a box OAR
            %
            % input:
            %   name:           Name of VOI as string
            %   dimensions:     Dimensions of the box as [x,y,z] array
            %
            % Name-Value pairs:
            %   'offset':       The offset of the VOI with respect to the
            %                   center of the geometry as [x,y,z] array
            %   'objectives':   Either a single objective or a cell array
            %                   of objectives
            %   'HU':           Houndsfield unit of the volume

            obj.volumes(end+1) = {matRad_PhantomVOIBox(name,'OAR',dimensions,varargin{:})};
            obj.updatecst();
        end

        function addSphericalOAR(obj,name,radius,varargin)
            % Adds a spherical OAR
            %
            % input:
            %   name:           Name of VOI as string
            %   radius:         Radius of the sphere
            %
            % Name-Value pairs:
            %   'offset':       The offset of the VOI with respect to the
            %                   center of the geometry as [x,y,z] array
            %   'objectives':   Either a single objective or a cell array
            %                   of objectives
            %   'HU':           Houndsfield unit of the volume

            obj.volumes(end+1) ={matRad_PhantomVOISphere(name,'OAR',radius,varargin{:})};
            obj.updatecst();
        end


        function [ct,cst] = getctcst(obj)
            %   Returns the ct and struct. The function also initializes
            %   the HUs in reverse order of defintion
            %
            % output
            %   ct:         matRad ct struct
            %   cst:        matRad cst struct
           
            %initialize the HU in reverse order of definition (objectives
            %defined at the start will have the highest priority in case of
            %overlaps)
            
            for i = 1:size(obj.cst,1)                    
                vIxVOI = obj.cst{end-i+1,4}{1};
                obj.ct.cubeHU{1}(vIxVOI) = 0; % assign HU 
            end
            
            ct  = obj.ct;
            cst = obj.cst;
        end

    end    

    methods (Access = private)

         function updatecst(obj)
            obj.cst =  obj.volumes{end}.initializeParameters(obj.ct,obj.cst);
        end

    end
end