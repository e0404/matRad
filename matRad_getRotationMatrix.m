function rotMat = matRad_getRotationMatrix(arg1,arg2,type,system)
% matRad function to return the rotation / transformation matrix for
% gantry and/or couch rotation. The Rotation matrix stands for a (1)
% counter-clockwise, (2) active rotation in the patient coordinate system
% that is performed on a (4) column vector (by premultiplying the matrix). 
% Per change of one of these directions a matrix transpose of the returned 
% matrix is required.
% 
% 
% call
%  rotMat = matRad_getRotationMatrix(gantryAngle,couchAngle,type,system)
%  rotMat = matRad_getRotationMatrix(pln,fieldIx,type,system)
%
% input
%   gantryAngle:    beam/gantry angle
%   couchAngle:     couch angle 
%   or  
%   pln:            matRad plan meta information struct
%   fieldIx:        index of the field / beam that is requested
%   
%   type:       optional parameter. Can be 'full' (both rotations as
%                   matrix, default) or 'gantry' / 'couch' for the 
%                   individual matrices
%   system:         optional coordinate system the transformation matrix is
%                   requested for. So far, only the default option 'LPS' is
%                   supported (right handed system).
%
% output
%   rotMat:         3x3 matrix that performs an active rotation around the 
%                   patient system origin via rotMat * x
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Parse arguments
%We need at least two and max 6 input arguments
narginchk(2,4);
%check if a pln and index have been given
if isstruct(arg1)
    pln = arg1;
    fieldIx = arg2;
    gantryAngle = pln.gantryAngles(fieldIx);
    couchAngle = pln.couchAngles(fieldIx);
else
    gantryAngle = arg1;
    couchAngle = arg2;
end 
% Coordinate System (only LPS so far)
if nargin < 4
    system = 'LPS';
end
%Gantry, Couch, or full transformation?
if nargin < 3
    type = 'full';
end

%% Set Up requested Rotation Matrix
switch system
    case 'LPS'
        %The LPS system is right-handed, gantry rotates counter-clockwise 
        %around z and couch rotates counter-clockwise around y
        
        %active, counter-clockwise rotation Matrix for Gantry around z 
        %with pre-multiplication of the matrix (R*x)
        %Note: Gantry rotation is physically an active rotation of a beam 
        %vector around the target / isocenterin the patient coordinate
        %system
        R_Gantry = [cosd(gantryAngle)  -sind(gantryAngle)  0; ...
            sind(gantryAngle)    cosd(gantryAngle)  0; ...
            0                            0           1];
        
        %active, counter-clockwise rotation for couch around y
        %with pre-multiplication of the matrix (R*x)
        %Note: Couch rotation is physically a passive rotation of the 
        %patient system around the beam target point / isocenter
        R_Couch = [cosd(couchAngle) 0 sind(couchAngle); ...
            0 1 0; ...
            -sind(couchAngle)    0    cosd(couchAngle)];
    otherwise
        error('matRad only supports LPS system so far');
end

switch type
    case 'full'
        %Since the couch rotation is physically a passive rotation of the
        %system, we perform it first and then rotate the gantry
        rotMat = R_Couch*R_Gantry;
    case 'gantry'
        rotMat = R_Gantry;
    case 'couch'
        rotMat = R_Couch;
    otherwise
        error(['Rotation matrix ''' type ''' not known, needs to be ''full'', ''gantry'' or ''couch''']);       
end

end

