function [ct,cst] =  matRad_createPhantomnew(varargin)
% matRad add margin function
%
% call
%   mVOIEnlarged = matRad_addMargin(mVOI,cst,vResolution,vMargin,bDiaElem)
%
% input
%   ctSize:         Array with size of the ct
%   ctResolution:   array with resolutions in mm
%   phantomType     Specify the phantom shape
%                   Possible values: 'cubic', 'sphere'
%   s
% output
%   ct:             matRad ct struct
%   cst:            matRad cst struct
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
% default values
defaultctSize = [200 200 50];
defaultResolutions = [3 3 3];
defaultCubicOAR = [100,100,25];%cubic target
defaultCubicPTV = [50,50,12.5];
defaultSphericalRadiusOAR = 40;%spherical target
defaultSphericalRadiusPTV = 12.5;
defaultCylindricalOAR =[40,25]; %cylindrical target
defaultCylindricalPTV = [12.5,12.5];
defaultShapeOAR = 'cubic';
defaultShapePTV = 'cubic';
expectedShapes = {'cubic','spheric','cylindric'};
%parse Input Variables
p = inputParser;
addOptional(p,'Resolutions',defaultResolutions);
addOptional(p,'ctSize', defaultctSize);
addOptional(p,'sphereRadiusOAR',defaultSphericalRadiusOAR);
addOptional(p,'sphereRadiusPTV',defaultSphericalRadiusPTV);
addOptional(p,'cubicOAR',defaultCubicOAR);
addOptional(p,'cubicPTV',defaultCubicPTV);
addOptional(p,'cylindricalOAR',defaultCylindricalOAR);
addOptional(p,'cylindricalPTV',defaultCylindricalPTV);
addParameter(p,'shapeOAR',defaultShapeOAR,...
             @(x) any(validatestring(x,expectedShapes)));
addParameter(p,'shapePTV',defaultShapePTV,...
             @(x) any(validatestring(x,expectedShapes)));
parse(p,varargin{:})
% CREATE CT-CUBE
ct.cubeDim      = p.Results.ctSize; % second cube dimension represents the x-coordinate
% set the resolution of x,y,z by looping over the structure
ct.resolution.x = p.Results.Resolutions(1);
ct.resolution.y = p.Results.Resolutions(2);
ct.resolution.z = p.Results.Resolutions(3);
ct.numOfCtScen  = 1;
% create a ct image series with zeros - it wil  l be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * -1000;
% CREATE VOI
ixOAR = 1;
ixPTV = 2;
% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';
cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';    
cst{ixPTV,3} = 'TARGET';
% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass  = 1;
cst{ixOAR,5}.alphaX       = 0.1000;
cst{ixOAR,5}.betaX        = 0.0500;
cst{ixOAR,5}.Priority     = 2;
cst{ixOAR,5}.Visible      = 1;
cst{ixOAR,5}.visibleColor = [1 0 0];
% define objective as struct for compatibility with GNU Octave I/O
cst{ixOAR,6}{1} = struct(DoseObjectives.matRad_SquaredOverdosing(10,30));
cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,5}.visibleColor = [0 0 1];
% define objective as struct for compatibility with GNU Octave I/O
cst{ixPTV,6}{1} = struct(DoseObjectives.matRad_SquaredDeviation(800,60));
% CREATE PHANTOM BASED OF CASE
%% Lets create either a cubic or a spheric phantom
% either 'cubic' or 'spheric'
% first the OAR
cubeHelper = zeros(ct.cubeDim);
switch p.Results.shapeOAR
   case {'cubic'}
      xLowOAR  = round(ct.cubeDim(2)/2 - p.Results.cubicOAR(1)/2);
      xHighOAR = round(ct.cubeDim(2)/2 + p.Results.cubicOAR(1)/2);
      yLowOAR  = round(ct.cubeDim(1)/2 - p.Results.cubicOAR(2)/2);
      yHighOAR = round(ct.cubeDim(1)/2 + p.Results.cubicOAR(2)/2);
      zLowOAR  = round(ct.cubeDim(3)/2 - p.Results.cubicOAR(3)/2);
      zHighOAR = round(ct.cubeDim(3)/2 + p.Results.cubicOAR(3)/2);
      display(zLowOAR)
      display(zHighOAR)
      for x = xLowOAR:1:xHighOAR
         for y = yLowOAR:1:yHighOAR
            for z = zLowOAR:1:zHighOAR
               cubeHelper(y,x,z) = 1;
            end
         end
      end
   case {'spheric'}
      radiusOAR = p.Results.sphereRadiusOAR;
      for x = 1:ct.cubeDim(2)
         for y = 1:ct.cubeDim(1)
            for z = 1:ct.cubeDim(3)
               currPost = [y x z] - round([ct.cubeDim./2]);
               if sqrt(sum(currPost.^2)) < radiusOAR
                  cubeHelper(y,x,z) = 1;
               end
            end
         end
      end
    case {'cylindric'} % cylinder aligned with z axis
      zLowOAR  = round(ct.cubeDim(3)/2 - p.Results.cylindricalOAR(1)/2);
      zHighOAR = round(ct.cubeDim(3)/2 + p.Results.cylindricalOAR(1)/2);
      display(zLowOAR)
      display(zHighOAR)
      cylindricalRadiusOAR = p.Results.cylindricalOAR(2)
      for x = 1:ct.cubeDim(2)
        for y = 1:ct.cubeDim(1)
            for z = zLowOAR:1:zHighOAR
                currPost = [y x z] - round([ct.cubeDim./2]);
                if sqrt(currPost(1)^2+currPost(2)^2) < cylindricalRadiusOAR
                   cubeHelper(y,x,z) = 1;
                end
           end
        end
      end
end
% extract the voxel indices and save it in the cst
cst{ixOAR,4}{1} = find(cubeHelper);
% second the PTV
cubeHelper = zeros(ct.cubeDim);
switch p.Results.shapePTV
   case {'cubic'}
      xLowPTV  = round(ct.cubeDim(2)/2 - p.Results.cubicPTV(1)/2);
      xHighPTV = round(ct.cubeDim(2)/2 + p.Results.cubicPTV(1)/2);
      yLowPTV  = round(ct.cubeDim(1)/2 - p.Results.cubicPTV(2)/2);
      yHighPTV = round(ct.cubeDim(1)/2 + p.Results.cubicPTV(2)/2);
      zLowPTV  = round(ct.cubeDim(3)/2 - p.Results.cubicPTV(3)/2);
      zHighPTV = round(ct.cubeDim(3)/2 + p.Results.cubicPTV(3)/2);
      cubeHelper = zeros(ct.cubeDim);
      for x = xLowPTV:1:xHighPTV
         for y = yLowPTV:1:yHighPTV
            for z = zLowPTV:1:zHighPTV
               cubeHelper(y,x,z) = 1;
            end
         end
      end
   case {'spheric'}
      radiusPTV = p.Results.sphereRadiusPTV;
      for x = 1:ct.cubeDim(2)
         for y = 1:ct.cubeDim(1)
            for z = 1:ct.cubeDim(3)
               currPost = [x y z] - round([ct.cubeDim./2]);
               if  sqrt(sum(currPost.^2)) < radiusPTV
                  cubeHelper(y,x,z) = 1;
               end
            end
         end
      end
  case {'cylindric'} % cylinder aligned with z axis
      zLowPTV = round(ct.cubeDim(3)/2 - p.Results.cylindricalPTV(1)/2);
      zHighPTV = round(ct.cubeDim(3)/2 + p.Results.cylindricalPTV(1)/2);
      cylindricalRadiusPTV = p.Results.cylindricalPTV(2);
      for x = 1:ct.cubeDim(2)
        for y = 1:ct.cubeDim(1)
            for z = zLowPTV:1:zHighPTV
                currPost = [y x z] - round([ct.cubeDim./2]);
                if sqrt(currPost(1)^2+currPost(2)^2) < cylindricalRadiusPTV
                   cubeHelper(y,x,z) = 1;
                end
           end
        end
      end
end
1
% extract the voxel indices and save it in the cst
cst{ixPTV,4}{1} = find(cubeHelper);
%Assign relative electron densities
vIxOAR = cst{ixOAR,4}{1};
vIxPTV = cst{ixPTV,4}{1};
ct.cubeHU{1}(vIxOAR) = 0; % assign HU of water
ct.cubeHU{1}(vIxPTV) = 0; % assign HU of water