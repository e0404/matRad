function matRad_MCNPBixelVisualization(ct, cst, stf, pln)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of bixel locations
%
% call
%   matRad_MCNPBixelVisualization(ct, cst, stf, pln)
%
% input
%       ct
%       cst
%       stf
%       pln
%
% output
%   ...
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 02/2019
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Selection of Gantry and Couch Angles

if pln.propStf.numOfBeams==1
    beamSelect = 1;     % index to spezify # of beam to plot
elseif pln.propStf.numOfBeams > 1
    for counterBeam = 1:pln.propStf.numOfBeams; beamList(counterBeam) = {['Gantry Angle: ', num2str(pln.propStf.gantryAngles(counterBeam)), ', Couch Angle: ', num2str(pln.propStf.couchAngles(counterBeam))]}; end
    [beamSelect,tf] = listdlg('ListString',beamList, 'Name','Please select a beam.', 'OKString', 'Select');
    if tf == 0 || length(beamSelect)~= 1
        warning('Please select only on combination of gantry and couch angle!')
        [beamSelect,tf] = listdlg('ListString',beamList, 'Name','Please select a beam.', 'OKString', 'Select');
    end
end

%% Plot geometry

% Prepare slice plot in isocenter plane
[x,y,z] = meshgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2), 1:ct.cubeDim(3));

% Get point on isocenter plane for given couch angle and plot isocenter
% plane through CT-data
if stf(beamSelect).couchAngle ~= 0 && stf(beamSelect).couchAngle ~= 90 && stf(beamSelect).couchAngle ~= 180 && stf(beamSelect).couchAngle ~= 270
    [xs, ys] = meshgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2));
    point = [round(stf(beamSelect).isoCenter(1)/ct.resolution.x) round(stf(beamSelect).isoCenter(2)/ct.resolution.y) round(stf(beamSelect).isoCenter(3)/ct.resolution.z)];
    normVec  = [sind(180 - (stf(beamSelect).couchAngle + 90)) 0 -cosd(180 - (stf(beamSelect).couchAngle + 90))];
    zs = ((xs-point(1))*normVec(1) + (ys-point(2))*normVec(2))/normVec(3) + point(3);
    
    slicePlot = slice(x, y, z, ct.cubeHU{1}, xs, ys, zs);
    %
    %     [xs, ys] = meshgrid(1:ct.cubeDim(1), 1:ct.cubeDim(2));
    %     point = [round(stf(beamSelect).isoCenter(2)/ct.resolution.y) round(stf(beamSelect).isoCenter(1)/ct.resolution.x) round(stf(beamSelect).isoCenter(3)/ct.resolution.z)];
    %     normVec  = [sind(180 - (stf(beamSelect).couchAngle + 90)) 0 -cosd(180 - (stf(beamSelect).couchAngle + 90))];
    %     zs = ((xs-point(1))*normVec(1) + (ys-point(2))*normVec(2))/normVec(3) + point(3);
    %
    %     slicePlot = slice(x, y, z, ct.cubeHU{1}, xs, ys, zs);
    
elseif stf(beamSelect).couchAngle == 0
    slicePlot = slice(x, y, z, ct.cubeHU{1}, [stf(beamSelect).isoCenter(1)/ct.resolution.x], [], []);
elseif stf(beamSelect).couchAngle == 90
    slicePlot = slice(x, y, z, ct.cubeHU{1}, [], [], [stf(beamSelect).isoCenter(3)/ct.resolution.z]);
elseif stf(beamSelect).couchAngle == 180
    slicePlot = slice(x, y, z, ct.cubeHU{1}, [], [stf(beamSelect).isoCenter(1)/ct.resolution.x], []);
elseif stf(beamSelect).couchAngle == 270
    slicePlot = slice(x, y, z, ct.cubeHU{1}, [], [], [stf(beamSelect).isoCenter(3)/ct.resolution.z]);
end


shading flat, colormap gray;
hold on

% Get target structures and plot isosurface
volumeDummy = zeros(ct.cubeDim);
for counterCst = 1:size(cst,1)
    if strcmp('TARGET', cst{counterCst,3})
        volumeDummy(cst{counterCst,4}{1}) = 1;
    end
end

isoPlot = patch(isosurface(x,y,z, volumeDummy, .9));

alpha(isoPlot, .7)
isoPlot.FaceColor = 'red';
isoPlot.EdgeColor = 'none';
%camlight('infinite')
lighting gouraud

% Plot source and bixel direction
for bixelCounter = 1:stf(beamSelect).numOfRays
    plot3( [stf(beamSelect).ray(bixelCounter).rayPosMLC(1)/ct.resolution.x (stf(beamSelect).isoCenter(1)+ stf(beamSelect).ray(bixelCounter).rayPos(1))/ct.resolution.x], ...
        [stf(beamSelect).ray(bixelCounter).rayPosMLC(2)/ct.resolution.y (stf(beamSelect).isoCenter(2)+ stf(beamSelect).ray(bixelCounter).rayPos(2))/ct.resolution.y], ...
        [stf(beamSelect).ray(bixelCounter).rayPosMLC(3)/ct.resolution.z (stf(beamSelect).isoCenter(3)+ stf(beamSelect).ray(bixelCounter).rayPos(3))/ct.resolution.z], ...
        '-', ...
        'Color', [135/256 206/256 250/256], ...
        'LineWidth',2.4)
end

% Plot MLC opening
for bixelCounter = 1:stf(beamSelect).numOfRays
    if stf(beamSelect).couchAngle ~= 0 && stf(beamSelect).couchAngle ~= 90 && stf(beamSelect).couchAngle ~= 180 && stf(beamSelect).couchAngle ~= 270
        
        xBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(1)-sind(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)-sind(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)+sind(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)+sind(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2];
        xBixel_MLC = [xBixel_MLC(1)/ct.resolution.x xBixel_MLC(2)/ct.resolution.y xBixel_MLC(3)/ct.resolution.z];
        yBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(2)-stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)+stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)-stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)+stf(beamSelect).bixelWidth/2];
        yBixel_MLC = [yBixel_MLC(1)/ct.resolution.x yBixel_MLC(2)/ct.resolution.y yBixel_MLC(3)/ct.resolution.z];
        zBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(3)-cosd(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)-cosd(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)+cosd(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)+cosd(stf(beamSelect).couchAngle)*stf(beamSelect).bixelWidth/2];
        zBixel_MLC = [zBixel_MLC(1)/ct.resolution.x zBixel_MLC(2)/ct.resolution.y zBixel_MLC(3)/ct.resolution.z];
    elseif stf(beamSelect).couchAngle == 0 || ...
            stf(beamSelect).couchAngle == 180 || ...
            stf(beamSelect).couchAngle == 90 || ...
            stf(beamSelect).couchAngle == 270
        
        xBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(1)...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(1)];
        xBixel_MLC = xBixel_MLC/ct.resolution.x;
        
        yBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(2)-stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)-stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)+stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(2)+stf(beamSelect).bixelWidth/2];
        yBixel_MLC = yBixel_MLC/ct.resolution.y;
        
        zBixel_MLC = [stf(beamSelect).ray(bixelCounter).rayPosMLC(3)-stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(3)+stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(3)+stf(beamSelect).bixelWidth/2 ...
            stf(beamSelect).ray(bixelCounter).rayPosMLC(3)-stf(beamSelect).bixelWidth/2];
        zBixel_MLC = zBixel_MLC/ct.resolution.z;
    end
    patch(xBixel_MLC, yBixel_MLC, zBixel_MLC, [135/256 206/256 250/256]); % 'b'); %[255/256 255/256 0/256])
end

axis tight
axis off


% Set rotation mode
cameratoolbar('SetMode','orbit')
cameratoolbar('SetCoordSys','y')
end