function stf = matRadBrachy_generateStf(ct,cst,pln)
%% config
matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating stf struct... ');
%% general stf information
stf = struct;
stf.radiationMode = 'brachy'; %pln.radiationMode
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

%translate to geometrix coordinates and save in stf

stf.targetVolume.Xvox = coordsX_vox; %hier fehlt noch eine Skalierung vo Vox auf realWorld-Koordinaten!
stf.targetVolume.Yvox = coordsY_vox;
stf.targetVolume.Zvox = coordsZ_vox;


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


%% generate 2D template shape
% Define steering file like struct.

[templXmesh,templYmesh] = meshgrid(1:8,1:6);
 
templX = reshape(templXmesh,[],1);
templY = reshape(templYmesh,[],1);
templZ = zeros(length(templX),1);
templ1 = ones(length(templX),1);

stf.template.template2D = [templX';templY';templZ';templ1'];
stf.template.Xscale = 5; %mm
stf.templateYscale = 5; %mm

%% oriention of template
%unit vectors of displaced, rotated template coordinate system
stf.orientation.Xdir = normalize([10,0,1],'norm');
stf.orientation.Ydir = normalize([0.02,1,-0.2],'norm');
stf.orientation.Zdir = cross(stf.orientation.Xdir,stf.orientation.Ydir);
stf.orientation.offset = [0,0,0]; %mm --> hier könnte das isocenter hin
%scale of X and Y direction in template
stf.orientation.Xscale = 50; %mm
stf.orientation.Yscale = 50; %mm

%throw error if directions are not orthogonal
assert(stf.orientation.Xdir*stf.orientation.Ydir' == 0,'Xdir and Ydir are not orthogonal')

% %% template rotation given by euler angles
% eulerPhi = 0;
% eulerTheta = 0;
% eulerPsi = 0;
% stf.template.rotAngles = [eulerPhi,eulerTheta,eulerPsi];
% stf.template.rotMatrix = eulerRot(eulerPhi,eulerTheta,eulerPsi);


%% needle
stf.needle.seedDistance = 20; %mm
stf.needle.seedsNo = 10;

%% generate seed positions
% seed positions can be generated from neeldes and template
% every seed will be indexed by (i, j) = (no. needle, no. seed in needle)

%preallocating space
stf.seeds.xPos = zeros(stf.needlesNo,stf.seedsPerNeedle);
stf.seeds.yPos = zeros(stf.needlesNo,stf.seedsPerNeedle);
stf.seeds.zPos = zeros(stf.needlesNo,stf.seedsPerNeedle);

%loop over all needles (i) and seeds (j)
for i = 1:stf.needlesNo
    for j = 1:stf.seedsPerNeedle
        
        xyzUntransformed = [stf.template.xPlanePositions(i);...
            stf.template.yPlanePositions(i); stf.needles(i,j)];
        xyzTransformed = stf.template.rotMatrix*xyzUntransformed;
        
        stf.seeds.xPos(i,j) = xyzTransformed(1);
        stf.seeds.yPos(i,j) = xyzTransformed(2);
        stf.seeds.zPos(i,j) = xyzTransformed(3);
    end
end

end

