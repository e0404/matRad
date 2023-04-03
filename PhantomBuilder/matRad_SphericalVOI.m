classdef  matRad_SphericalVOI < matRad_VOIVolume   
    % matRad_SphericalVOI implements a class that helps to create spheric VOIs
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
    properties
        radius;
    end
    methods
        function obj = matRad_SphericalVOI(name,type,radius,varargin)
            p = inputParser;
            addOptional(p,'objectives',{});
            addOptional(p,'offset',[0,0,0]);
            parse(p,varargin{:});

            obj@matRad_VOIVolume(name,type,p); %call superclass constructor
            obj.radius = radius;
        end
        
        function [cst] = initializeParameters(obj,ct,cst)
            cst = initializeParameters@matRad_VOIVolume(obj,cst);
            center = round([ct.cubeDim/2]);
            VOIHelper = zeros(ct.cubeDim);
            offsets = obj.offset;

            for x = 1:ct.cubeDim(2)
                for y = 1:ct.cubeDim(1)
                   for z = 1:ct.cubeDim(3)
                      currPost = [x y z]  + offsets - center;
                      if  (sqrt(sum(currPost.^2)) < obj.radius)
                            VOIHelper(y,x,z) = 1;
                      end
                   end
                end
            end
            
            cst{obj.idx,4}{1} = find(VOIHelper);
            
        end
    end
end