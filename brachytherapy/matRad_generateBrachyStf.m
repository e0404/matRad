function stf = matRad_generateBrachyStf(ct,cst,pln)
% matRad_generateBrachyStf generates matRad steering information generation for brachy
% 
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct

%% config
matRad_cfg = MatRad_Config.instance();
addpath(fullfile( matRad_cfg.matRadRoot));
matRad_cfg.dispInfo('matRad: Generating stf struct... ');

if ~isfield(pln,'propStf')
matRad_cfg.dispError('no applicator information in pln struct');
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

%translate to geometric coordinates and save in stf

stf.targetVolume.Xvox = ct.x(coordsX_vox); % angabe in mm
stf.targetVolume.Yvox = ct.y(coordsY_vox);
stf.targetVolume.Zvox = ct.z(coordsZ_vox);


% % calculate rED or rSP from HU
% ct = matRad_calcWaterEqD(ct, pln);

% take only voxels inside patient
V = [cst{:,4}];
V = unique(vertcat(V{:}));

% ignore densities outside of contours
eraseCtDensMask = ones(prod(ct.cubeDim),1);
eraseCtDensMask(V) = 0;
for i = 1:ct.numOfCtScen
    ct.cube{i}(eraseCtDensMask == 1) = 0;
end

        %% save VOI coordinates
        %translate to geometric coordinates and save in stf
        stf.targetVolume.Xvox = ct.x(coordsX_vox);
        stf.targetVolume.Yvox = ct.y(coordsY_vox);
        stf.targetVolume.Zvox = ct.z(coordsZ_vox);
%         stf.targetVolume.Xvox = ct.x(coordsX_vox); % given in mm
%         stf.targetVolume.Yvox = ct.y(coordsY_vox);
%         stf.targetVolume.Zvox = ct.z(coordsZ_vox);
        %% copy meta info from pln
        stf.radiationMode = pln.radiationMode;
        stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
        stf.numOfNeedles = pln.propStf.template.numOfXPoints*pln.propStf.template.numOfYPoints; %might be changed for coordinate based pln input
        stf.totalNumOfBixels = stf.numOfSeedsPerNeedle*stf.numOfNeedles; % means total number of seeds 

        %% generate 2D template shape
        % in case, the template was defined as grid, 4xN position matrix is
        % produced. Each column is the 3+1 vector of one template hole (fourth entry always one, z coord - third entry zero)
        % the template origin is set at its center.
        if strcmp(pln.propStf.templateDirect, 'none')
            xNum = pln.propStf.template.numOfXPoints;
            yNum = pln.propStf.template.numOfYPoints;
            xSc = pln.propStf.template.xScale;
            ySc = pln.propStf.template.yScale;
            [templXmesh,templYmesh] = meshgrid(-(xNum-1)/2:(xNum-1)/2,-(yNum-1)/2:(yNum-1)/2);
            templX = xSc*reshape(templXmesh,[],1);
            templY = ySc*reshape(templYmesh,[],1);
            templZ = zeros(length(templX),1);
            templ1 = ones(length(templX),1);
            stf.template.template2D = [templX';templY';templZ';templ1'];
        else
            stf.template.template2D = pln.propStf.templateDirect;
        end
        
        %% generate seed positions
        % seed positions can be generated from neeldes, template and oriantation
        % needles are assumed to go trough the template vertically
        
        % transformed template
        stf.template.template3D = pln.propStf.shiftRotMtx*stf.template.template2D;

        % needle position
        d = pln.propStf.needle.seedDistance;
        seedsNo = pln.propStf.needle.seedsNo;
        needleDist(1,1,:) = d.*[1:seedsNo]'; % 1x1xN Array with seed positions on needle
        needleDir = needleDist.*pln.propStf.shiftRotMtx(1:3,3);
        seedPos_coord_need_seed = needleDir + stf.template.template3D(1:3,:);
        seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed,1);
        % the output array has the dimentions (needleNo,seedNo,coordinates)
        X = seedPos_need_seed_coord(:,:,1);
        Y = seedPos_need_seed_coord(:,:,2);
        Z = seedPos_need_seed_coord(:,:,3);
        
        stf.seedPoints.x = reshape(X,1,[]);
        stf.seedPoints.y = reshape(Y,1,[]);
        stf.seedPoints.z = reshape(Z,1,[]);
        
        matRad_cfg.dispInfo('...100% ');

        % trow warning if seed points are more then twice the central
        % distange outsidethe TARGET volume or if no sed points are in the
        % target volume
        
        if (max(stf.seedPoints.x) >= 2*max(stf.targetVolume.Xvox) || min(stf.seedPoints.x) <= 2*min(stf.targetVolume.Xvox) ||...
            max(stf.seedPoints.y) >= 2*max(stf.targetVolume.Yvox) || min(stf.seedPoints.y) <= 2*min(stf.targetVolume.Yvox) || ...
            max(stf.seedPoints.z) >= 2*max(stf.targetVolume.Zvox) || min(stf.seedPoints.z) <= 2*min(stf.targetVolume.Zvox))
                matRad_cfg.dispWarning('Seeds far outside the target volume');
        end
        if (max(stf.targetVolume.Xvox) <= min(stf.seedPoints.x) || min(stf.targetVolume.Xvox) >= max(stf.seedPoints.x) ||...
            max(stf.targetVolume.Yvox) <= min(stf.seedPoints.y) || min(stf.targetVolume.Yvox) >= max(stf.seedPoints.y) ||...
            max(stf.targetVolume.Zvox) <= min(stf.seedPoints.z) || min(stf.targetVolume.Zvox) >= max(stf.seedPoints.z))
                matRad_cfg.dispWarning('no seed points in VOI')
        end    
end

