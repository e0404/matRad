function [dosePoints,volume,patch] = generateNormalTissue(obj,rpItem,additionalInfo)
%GENERATENORMALTISSUE Summary of this function goes here
%   Detailed explanation goes here
global cm mm


if nargin<3
    additionalInfo = [];
end
if ~isfield(additionalInfo,'margin')
    additionalInfo.margin = 0.75*cm;
end
if ~isfield(additionalInfo,'samplingDensity')
    additionalInfo.samplingDensity = 1/mm^3;
end

if ~isfield(additionalInfo,'samplingDensity')
    additionalInfo.maxPoints = 1500;
end

dwellPositions = cell2mat(obj.RP.(rpItem).seeds.dwellPosition);


rsNames        = fieldnames(obj.RS);

ir = 0;
found = 0;
while ir<length(rsNames) && ~found
    ir    = ir+1;
    found = strcmp(obj.RS.(rsNames{ir}).additionalInfo.structID,obj.RP.(rpItem).additionalInfo.structID);
end

item = rsNames{ir};
rois = obj.RS.(item).structures;

fnames = fieldnames(rois);
center = mean(dwellPositions,1);

%dwellPositions-center*
dist = dwellPositions-ones(size(dwellPositions,1),1)*center;
dist = max(sqrt(sum(dist.*dist,2)));


%volume        = 4/3*pi*(dist+additionalInfo.margin)^3;

additionalInfo.margin        = dist+additionalInfo.margin;
minXYZ        = center-(additionalInfo.margin);
maxXYZ        = center+(additionalInfo.margin);

volume        = 4/3*pi*(additionalInfo.margin)^3;
noPointsInVol = round(volume*additionalInfo.samplingDensity);


% diffVol = volume;
% for in=1:length(fnames)
%     diffVol   = diffVol-rois.(fnames{in}).volume;
% end


if isfield(additionalInfo,'maxPoints')
    if additionalInfo.maxPoints<noPointsInVol
        noPointsInVol = 2*additionalInfo.maxPoints;
    end
end

% generate random dose Points
points = sobolset(3);
points = scramble(points,'MatousekAffineOwen');
points = points(1:noPointsInVol,:);
points  = ones(noPointsInVol,1)*minXYZ+(ones(noPointsInVol,1)*(maxXYZ-minXYZ)).*points;

dist = points-ones(noPointsInVol,1)*center;
dist = sqrt(sum(dist.*dist,2));

points = points(dist<=additionalInfo.margin,:);

if isfield(additionalInfo,'coordinates')
    
    indAcc = min(additionalInfo.coordinates.x)<=points(:,1) & points(:,1)<=max(additionalInfo.coordinates.x) &...
        min(additionalInfo.coordinates.y)<=points(:,2) & points(:,2)<=max(additionalInfo.coordinates.y) &...
        min(additionalInfo.coordinates.z)<=points(:,3) & points(:,3)<=max(additionalInfo.coordinates.z);
    points = points(indAcc,:);
    
end

ind = false(size(points,1),length(fnames));
h = waitbar(0,'Please wait generating normal tissue...');
for in=1:length(fnames)  
    tempInd = false(size(points,1),length(rois.(fnames{in}).outerContourMap));
    for im=1: length(rois.(fnames{in}).outerContourMap)
        [~,tempInd(:,im)] = getPointsInVolume(rois.(fnames{in}).outerContourMap(im),points);
        
    end
    ind(:,in) = any(tempInd,2);
    waitbar(in / length(fnames))
end
close(h)

ind = prod(~ind,2)==1;

dosePoints = points(ind,:);

if nargout>2

[X,Y,Z] = sphere;
fvc = surf2patch(X,Y,Z,'triangles');
patch.faces = fvc.faces;
patch.vertices = fvc.vertices*additionalInfo.margin+ones(size(fvc.vertices,1),1)*mean(dwellPositions,1);
end

end

