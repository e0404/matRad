%% Euler Rotation
% Given the euler angles, this function generates the euler rotation
% matrix.
% Explanation of euler angles: https://en.wikipedia.org/wiki/Euler_angles
% In short, the system is turned around the z then x' then z'' axis.

% Positive angles -> turn axis
% Negative angles -> coordinate transformation


function RotMtx = eulerRot(phi,theta,psi)

if nargin == 1 % in that case, the angles are probably given as vector
    theta = phi(2);
    psi = phi(3);
    phi = phi(1);
end

% Rotation matrces
R1 = [cos(phi),-sin(phi),0 ; sin(phi),cos(phi),0 ; 0,0,1];
R2 = [1,0,0 ; 0,cos(theta),-sin(theta) ; 0,sin(theta),cos(theta)];
R3 = [cos(psi),-sin(psi),0 ; sin(psi),cos(psi),0 ; 0,0,1];

% Rotation
RotMtx = R1*R2*R3;


end
