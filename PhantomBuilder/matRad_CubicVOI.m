classdef matRad_CubicVOI < matRad_VOIVolume
    % matRad_CubicVOI implements a class that helps to create cubic VOIs
    %
    % References 
    %     -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    properties %additional property of cubic objects 
        dimensions;
    end
    methods
        function obj = matRad_CubicVOI(name,type,dimensions,varargin)
            p = inputParser; 
            addOptional(p,'objectives',{});
            addOptional(p,'offset',[0,0,0]);
            parse(p,varargin{:});

            obj@matRad_VOIVolume(name,type,p); %call superclass constructor
            obj.dimensions = dimensions;
        end
        
        function [cst] = initializeParameters(obj,ct,cst)
            cst = initializeParameters@matRad_VOIVolume(obj,cst);
            center = round(ct.cubeDim/2);
            VOIHelper = zeros(ct.cubeDim);
            offsets = obj.offset;
            dims = obj.dimensions;

            for x = center(1)+offsets(1) - dims(1)/2 :1: center(1) + offsets(1) + dims(1)/2
                for y = center(2)+offsets(2) - dims(2)/2 :1: center(2) + offsets(2) + dims(2)/2
                   for z = round(center(3)+offsets(3) - dims(3)/2) :1: center(3) + offsets(3) + dims(3)/2
                        VOIHelper(y,x,z) = 1;
                        
                   end
                end
            end
            
            cst{obj.idx,4}{1} = find(VOIHelper);
            

        end
    end
end  