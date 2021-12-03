function stf = matRad_generateBrachyStf(ct,cst,pln, visMode)
% matRad_generateBrachyStf generates matRad steering information generation for brachy
%   will be called within matRad_generateStf if radiation mode is 'brachy'
%
% call
%   stf = matRad_generateBrachyStf(ct,cst,pln,visMode) 
%   
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this 
%               value to 1,2,3 (optional)
%
% output
%   stf:        matRad steering information struct
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%% meta info from pln
stf.radiationMode = pln.radiationMode;
stf.numOfSeedsPerNeedle = pln.propStf.needle.seedsNo;
stf.numOfNeedles = nnz(pln.propStf.template.activeNeedles);
stf.totalNumOfBixels = stf.numOfSeedsPerNeedle*stf.numOfNeedles; % means total number of seeds 

%% generate 2D template points
% the template origin is set at its center. In the image coordinate system,
% the center will be positioned at the bottom of the volume of interest.

[row,col] = find(pln.propStf.template.activeNeedles);
templX = col*pln.propStf.bixelWidth + pln.propStf.templateRoot(1) - (13+1)/2*pln.propStf.bixelWidth;
templY = row*pln.propStf.bixelWidth + pln.propStf.templateRoot(2) - (13+1)/2*pln.propStf.bixelWidth;
templZ = ones(size(col))                 + pln.propStf.templateRoot(3);


stf.template = [templX';templY';templZ'];



        
%% generate seed positions
% seed positions can be generated from neeldes, template and oriantation
% needles are assumed to go trough the template vertically

% needle position
d = pln.propStf.needle.seedDistance;
seedsNo = pln.propStf.needle.seedsNo;
needleDist(1,1,:) = d.*[0:seedsNo-1]'; % 1x1xN Array with seed positions on needle
needleDir = needleDist.*[0;0;1];
seedPos_coord_need_seed = needleDir + stf.template;
seedPos_need_seed_coord = shiftdim(seedPos_coord_need_seed,1);
% the output array has the dimentions (needleNo,seedNo,coordinates)
X = seedPos_need_seed_coord(:,:,1);
Y = seedPos_need_seed_coord(:,:,2);
Z = seedPos_need_seed_coord(:,:,3);

stf.seedPoints.x = reshape(X,1,[]);
stf.seedPoints.y = reshape(Y,1,[]);
stf.seedPoints.z = reshape(Z,1,[]);

matRad_cfg.dispInfo('...100% ');

%%visualize results of visMode is nonzero
% plot 3D seed positions
if visMode > 0
    clf
    SeedPoints = plot3(stf.seedPoints.x,stf.seedPoints.y,stf.seedPoints.z,'.','DisplayName', 'seed points','Color','black','markersize',5);
    title( '3D Visualization of seed points')
    xlabel('X (left) [mm]')
    ylabel('Y (posterior) [mm]')
    zlabel('Z (superior) [mm]')
    hold on
    

    % plot 3d VOI points
    TargX = stf.targetVolume.Xvox;
    TargY = stf.targetVolume.Yvox;
    TargZ = stf.targetVolume.Zvox;
    %Prostate = plot3(TargX,TargY,TargZ,'.', 'Color','b','DisplayName', 'prostate');
    
    P = [TargX',TargY',TargZ'];
    k = boundary(P,1);
    trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','red','FaceAlpha',0.1,'LineStyle','none')
    hold off;
end


% trow warning if seed points are more then twice the central
% distange outsidethe TARGET volume or if no sed points are in the
% target volume

if (max(stf.seedPoints.x-pln.propStf.templateRoot(1)) >= 4*max(stf.targetVolume.Xvox-pln.propStf.templateRoot(1)) ||...
        min(stf.seedPoints.x-pln.propStf.templateRoot(1)) <= 4*min(stf.targetVolume.Xvox-pln.propStf.templateRoot(1)) ||...
    max(stf.seedPoints.y-pln.propStf.templateRoot(2)) >= 4*max(stf.targetVolume.Yvox-pln.propStf.templateRoot(2)) ||...
        min(stf.seedPoints.y-pln.propStf.templateRoot(2)) <= 4*min(stf.targetVolume.Yvox-pln.propStf.templateRoot(2)) || ...
    max(stf.seedPoints.z-pln.propStf.templateRoot(3)) >= 4*max(stf.targetVolume.Zvox-pln.propStf.templateRoot(3)) ||...
        min(stf.seedPoints.z-pln.propStf.templateRoot(3)) <= 4*min(stf.targetVolume.Zvox-pln.propStf.templateRoot(3)))
        matRad_cfg.dispWarning('Seeds far outside the target volume');
end
if (max(stf.targetVolume.Xvox) <= min(stf.seedPoints.x) || min(stf.targetVolume.Xvox) >= max(stf.seedPoints.x) ||...
    max(stf.targetVolume.Yvox) <= min(stf.seedPoints.y) || min(stf.targetVolume.Yvox) >= max(stf.seedPoints.y) ||...
    max(stf.targetVolume.Zvox) <= min(stf.seedPoints.z) || min(stf.targetVolume.Zvox) >= max(stf.seedPoints.z))
        matRad_cfg.dispWarning('no seed points in VOI')
end    

end

