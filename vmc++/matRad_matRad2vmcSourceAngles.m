function angles = matRad_matRad2vmcSourceAngles(gantryAngle,couchAngle)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad convert gantry and couch angles to angles used by vmc++
%
% call
%   matRad_matRad2mvcSourceAngles(gantryAngle,couchAngle)
%
% input
%   gantryAngle:    gantry angle
%   couchAngle:     couch angle
%
%
% References
%   Notes 25 July 2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch -cosd(couchAngle).*sind(gantryAngle)
    % do special cases first, when cosd(thetaY) == 0
    % -cosd(couchAngle)*sind(gantryAngle) = sind(thetaY)
    
    case 1
        thetaY = 90;
        thetaZ = 0;
        
        if couchAngle == 0
            % then gantryAngle == 270
            thetaX = 90;
        else
            % then couchAngle == 180, gantryAngle == 90
            thetaX = 270;
        end
        
    case -1
        thetaY = 270;
        thetaZ = 0;
        if couchAngle == 0
            % then gantryAngle == 90
            thetaX = 90;
        else
            % then couchAngle == 180, gantryAngle == 270
            thetaX = 270;
        end
        
    otherwise
        % general case, cosd(thetaY) ~= 0
        
        % first determine thetaY; note that we may have to take
        % supplementary angle later
        thetaY = asind(-cosd(couchAngle).*sind(gantryAngle));
        
        % now determine thetaX and thetaZ from the x and y components
        thetaX = atan2d(cosd(gantryAngle)./cosd(thetaY),sind(couchAngle).*sind(gantryAngle)./cosd(thetaY));
        thetaZ = atan2d(-sind(couchAngle)./cosd(thetaY),cosd(couchAngle).*cosd(gantryAngle)./cosd(thetaY));
        
        % verify that the remaining angular relations are satisfied
        % if not, take thetaY = 180-thetaY and recalculate thetaX and
        % thetaZ
        if ~verifiedRelations(gantryAngle,couchAngle,thetaX,thetaY,thetaZ)
            thetaY = 180-thetaY;
            thetaX = atan2d(cosd(gantryAngle)./cosd(thetaY),sind(couchAngle).*sind(gantryAngle)./cosd(thetaY));
            thetaZ = atan2d(-sind(couchAngle)./cosd(thetaY),cosd(couchAngle).*cosd(gantryAngle)./cosd(thetaY));
        end
        
end

% now verify for a final time that the remaining angular relations are satisfied
if ~verifiedRelations(gantryAngle,couchAngle,thetaX,thetaY,thetaZ)
    error('Angular relations are not satisfied for some reason');
end

angles = [thetaX thetaY thetaZ];

end


function ver = verifiedRelations(gantryAngle,couchAngle,thetaX,thetaY,thetaZ)

ver = true;

% need to verify four relations
% if any of them fail, then fail the test
if abs(sind(gantryAngle) + sind(thetaX).*sind(thetaY).*cosd(thetaZ) + cosd(thetaX).*sind(thetaZ)) > eps*10e3
    ver = false;
end

if abs(-sind(couchAngle).*cosd(gantryAngle) + cosd(thetaX).*sind(thetaY).*cosd(thetaZ) - sind(thetaX).*sind(thetaZ)) > eps*10e3
    ver = false;
end

if abs(sind(thetaX).*sind(thetaY).*sind(thetaZ)-cosd(thetaX).*cosd(thetaZ)) > eps*10e3
    ver = false;
end

if abs(-cosd(couchAngle)+cosd(thetaX).*sind(thetaY).*sind(thetaZ)+sind(thetaX).*cosd(thetaZ)) > eps*10e3
    ver = false;
end

end