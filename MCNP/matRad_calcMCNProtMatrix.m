function rotMatrix = matRad_calcMCNProtMatrix(gantryAngle, couchAngle)
% Calculate MCNP rotation matrix
% 
% call
%   rotMatrix = matRad_calcMCNProtMatrix(gantryAngle, couchAngle)
%
% input
%   gantryAngle
%   couchAngle
%
% output
%   rotMatrix:      Roation matrix as MCNP input.
%                   Note: Rotation angle for couch rotation in opposite
%                   direction for matRad than for usual rotation in eD.
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Phi = gantryAngle;
Theta = couchAngle;

matrixGantryRot_z =  @(Phi)     [cosd(Phi) -sind(Phi) 0; ...
                                sind(Phi) cosd(Phi) 0; ...
                                0 0 1];

matrixCouchRot_y = @(Theta)     [cosd(360-Theta) 0 sind(360-Theta);...
                                0 1 0; ...
                                -sind(360-Theta) 0 cosd(360-Theta)];

getCosAngle = @(Phi, Theta, vec1, vec2) ...
    vec1*(matrixCouchRot_y(Theta)*matrixGantryRot_z(Phi)*vec2')/...
    (sqrt(sum([vec1].^2))*sqrt(sum([(matrixCouchRot_y(Theta)*matrixGantryRot_z(Phi)*vec2')].^2)));

e_x = [1 0 0];
e_y = [0 1 0];
e_z = [0 0 1];

rotMatrix = [getCosAngle(Phi, Theta, e_x, e_x) getCosAngle(Phi, Theta,e_y , e_x) getCosAngle(Phi, Theta, e_z, e_x);...
            getCosAngle(Phi, Theta, e_x, e_y) getCosAngle(Phi, Theta, e_y, e_y) getCosAngle(Phi, Theta, e_z, e_y); ...
            getCosAngle(Phi, Theta, e_x, e_z) getCosAngle(Phi, Theta, e_y, e_z) getCosAngle(Phi, Theta, e_z, e_z)];