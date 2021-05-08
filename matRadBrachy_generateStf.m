function stf = matRadBrachy_generateStf(ct,cst,pln,visMode)
% matRad Brachytherapy steering information generation
% This function generates:
%       The Volume of interest 
% 
% call
%   stf = matRadBrachy_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
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

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating stf struct... ');

if nargin < 4
    visMode = 0;
end

%% generate image coordinates

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;vertcat(cst{i,4}{:})];
    end
end

% Remove double voxels
V = unique(V);
% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;
    
% add margin
addmarginBool = matRad_cfg.propStf.defaultAddMargin;
if isfield(pln,'propStf') && isfield(pln.propStf,'addMargin')
   addmarginBool = pln.propStf.addMargin; 
end

if addmarginBool
    voiTarget = matRad_addMargin(voiTarget,cst,ct.resolution,ct.resolution,true);
    V   = find(voiTarget>0);
end

% throw error message if no target is found
if isempty(V)
    matRad_cfg.dispError('Could not find target.');
end

% Convert linear indices to 3D voxel coordinates
[coordsY_vox, coordsX_vox, coordsZ_vox] = ind2sub(ct.cubeDim,V);

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep 'basedata' filesep fileName]);
   SAD = machine.meta.SAD;
catch
   matRad_cfg.dispError('Could not find the following machine file: %s',fileName); 
end

if strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
      
    availableEnergies = [machine.data.energy];
    availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
    
    if sum(availablePeakPos<0)>0
       matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file') 
    end
    %clear machine;
end

% calculate rED or rSP from HU
ct = matRad_calcWaterEqD(ct, pln);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

%% generate needle information
% Define steering file like struct. Prellocating for speed.
stf = struct;

% position of points
xPointsNo = 20;
yPointsNo = 20;
pointsNo = xPointsNo*yPointsNo;
needlesNo = pointsNo;
pointDistance = 10; %mm %kp, ob das stimmt...

for i = 1:pointsNo
stf.template.xPlanePositions(i,1) = pointDistance*(mod(i-1,xPointsNo)+1);
stf.template.yPlanePositions(i,1) = pointDistance*(ceil(xPointsNo*i/pointsNo));
end

% template offset(position)
stf.template.offset = [30;30;30]; %mm

% template rotation 
%given by euler angles
eulerPhi = pi/4;
eulerTheta = pi/4;
eulerPsi = pi/4;
stf.template.rotAngles = [eulerPhi,eulerTheta,eulerPsi];
stf.template.rotMatrix = eulerRot(eulerPhi,eulerTheta,eulerPsi);


%% needle seeds
%Diameter
stf.needleDiameter = 0.01; %mm %not used for first prototype, point sources assumed

%in the following matrix, lines represent needle indices, colums seed
%indices.
seedDistance = 1; %mm
seedsNo = 10;

for i = 1:pointsNo, j = 1:seedsNo;
    stf.needles(i,j) = seedDistance*j;
end

%% generate seed positions
% seed positions can be generated from neeldes and template
% every seed will be indexed by (i, j) = (no. needle, no. seed in needle)

%preallocating space
stf.seeds.xPos = zeros(needlesNo,seedsNo);
stf.seeds.xPos = zeros(needlesNo,seedsNo);
stf.seeds.xPos = zeros(needlesNo,seedsNo);

%loop over all needles (i) and seeds (j)
for i = 1:neeldeNo
    for j = 1:seedsNo
        
        xyzUntransformed = [stf.template.xPlanePositions(i);...
            stf.template.xPlanePositions(i); stf.needles(i,j)];
        xyzTransformed = stf.template.rotMatrix*xyzUntransformed;
        
        stf.seeds.xPos(i,j) = xyzTransformed(1);
        stf.seeds.yPos(i,j) = xyzTransformed(2);
        stf.seeds.zPos(i,j) = xyzTransformed(3);
    end
end

end

