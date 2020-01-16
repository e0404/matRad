%%%% script for preparing box phantoms with slabs - prepared for Paul

clear, clc

% settings

% water box phantom size and resolution
vars.boxSize = 160 * ones(1,3);
vars.res = 2 * ones(1,3);

vars.geoShape = "Rectangle"; 

% size of the slab and alignment

slabYs = 3;     %0:3;   % size of slab in y axis
slabXs = 15;    %11:18; % use this to take the slab in and out of the beam
                        % as we discussed, there is no size parameter in x axis
slabZs = 3;     %0:3;   % size of slab in z axis


alignmentsX =  -35;     %-35:-10; % alignment in the beam direction(slab up and down)
                       % -10 so that it does not stand in the bragg peak
                       % lower range depends on the energy of the particle
                       

vars.alignment = [0 -15 0]; % hard coded stuff!
vars.alignment(1) = alignmentsX;

% as I said, not very smart way to approach the problem!

% vars.geoSize(1) = randsample(slabYs,1); 
% vars.geoSize(2) = randsample(slabXs,1); 
% vars.geoSize(3) = randsample(slabZs,1);
vars.geoSize = [slabYs slabXs slabZs];


vars.slab_sp = 2; % stopping power of the slab
vars.tissue_sp = 1.1; % stopping power of the tissue

[ct, cst] = makeBoxphantom(vars.boxSize, vars.res, vars.tissue_sp);

mask_zeros = zeros(ct.cubeDim);
mask_zeros(ct.cube{1} == 0) = 1;

pln.radiationMode = 'protons';
pln.machine = 'generic_TOPAS_cropped';

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = 0;
pln.propStf.couchAngles   = 0;
pln.propStf.bixelWidth    = 150;
pln.propStf.longitudinalSpotSpacing = 1.5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;


pln.propOpt.bioOptimization = 'none';

% reading the isocenter
isoCenter = pln.propStf.isoCenter; 

isoCenter(1) = floor(isoCenter(1)/ct.resolution.y);
isoCenter(2) = floor(isoCenter(2)/ct.resolution.x);
isoCenter(3) = floor(isoCenter(3)/ct.resolution.z);
% turning it into voxels and aligning it

vars.slab_loc = isoCenter + vars.alignment;

% bulding the mask for where the slab is
mask = slabGeometry(vars, ct.cubeDim);

ct.cube{1}(mask == 1) = vars.slab_sp;
ct.cube{1}(mask_zeros == 1) = 0;
ct = matRad_electronDensitiesToHU(ct);

contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),3,'color','black');
