function XMLplan = matRad_create_XMLplan(pln,stf,w,MCsettings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create structure to export beam sport for MC simulation
% 
% call
%   XMLplan = matRad_create_XMLplan()
% 
% input
%   pln:            matRad plan
%   stf:            matRad steering information
%   w (optional):   precalculated weights
%
% output
%   XMLplan:                plan structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('generating XMLplan structure')

if ~isfield(MCsettings,'minRelWeight')
    MCsettings.minRelWeight = 0.0;
end

maxParticlesInSpot = 1e6*max(w(:));
minParticlesInSpot = 1e6*min(w(:));
minNbParticlesSpot = round(max([MCsettings.minRelWeight*maxParticlesInSpot,1]));

% Compute bixel index as it is done on matRad_calParticleDose
counter = 0;

numOfBeams = [];
if isfield(pln,'propStf')
  numOfBeams = pln.propStf.numOfBeams;
else
  numOfBeams = pln.numOfBeams;
end

for i = 1:numOfBeams % loop over all beams

    for j = 1:stf(i).numOfRays % loop over all rays

        if ~isempty(stf(i).ray(j).energy)

            for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                counter = counter + 1;

                % remember beam and  bixel number
                tmp.beamNum(counter)  = i;
                tmp.rayNum(counter)   = j;
                tmp.bixelNum(counter) = k;
            end
        end
    end
end

% load machine file, required to read machine.data.initFocus.SisFWHMAtIso
fileName = [pln.radiationMode '_' pln.machine];
try
   load(fileName);
catch
   error(['Could not find the following machine file: ' fileName ]);
end
availableEnergies = [machine.data.energy];

for beamNb = 1:numOfBeams % loop over beams

  XMLplan(beamNb).nbParticles = 0;
  XMLplan(beamNb).beamNb = beamNb;
  XMLplan(beamNb).couchAngle = stf(beamNb).couchAngle;
  XMLplan(beamNb).gantryAngle = stf(beamNb).gantryAngle;
  if strcmp(pln.radiationMode,'protons')
    XMLplan(beamNb).ion = '1H';
    XMLplan(beamNb).ionZ=1;
    XMLplan(beamNb).ionA=1;

  elseif strcmp(pln.radiationMode,'carbon')
    XMLplan(beamNb).ion = '12C';
    XMLplan(beamNb).ionZ=6;
    XMLplan(beamNb).ionA=12;
  end

  % helper function for energy selection
  round2 = @(a,b)round(a*10^b)/10^b;

  % Find IES values
  iesArray = [];
  for rayNb = 1:stf(beamNb).numOfRays
    iesArray = unique([iesArray stf(beamNb).ray(rayNb).energy]);
  end

  iesNb = 0; % used to get rid off IES in the iesArray with no particles (weigth = 0)
  for iesArrayIx = 1:length(iesArray) % find focus and "Voxel (x,y,nbParticles)" for each IES
    clear iesFocus iesEnergy;
    newIES = true;
    iesEnergy = iesArray(iesArrayIx);
    iesFocus = 0;

    for rayNb = 1:stf(beamNb).numOfRays % loop over rays

      % find index of used energy (round to keV for numerical reasons
      bixelNb = find(round2(stf(beamNb).ray(rayNb).energy,4) == round2(iesEnergy,4)==1);
      energyIx = find(round2(iesEnergy,4) == round2(availableEnergies,4)==1);

      if length(bixelNb)==1 % one IES found

        bixelIndex = find([tmp.beamNum==beamNb & tmp.rayNum==rayNb & tmp.bixelNum==bixelNb]==1);

        voxel_nbParticles = w(bixelIndex);
        voxel_nbParticles = round(1e6*voxel_nbParticles);

        % check whether there are (enough) particles for beam delivery
        if (voxel_nbParticles>minNbParticlesSpot)

          rayPos_bev = stf(beamNb).ray(rayNb).rayPos_bev;
          % matRad, x-y(beam)-z
          % HITXML, x-y-z(beam)
          %
          %        -X
          %        |
          %  -Y -- o -- +Y
          %        |
          %        +X
          %
          voxel_x = -rayPos_bev(3);
          voxel_y = rayPos_bev(1);

          rayIESenergy=stf(beamNb).ray(rayNb).energy(bixelNb);
          rayIESfocusIx=stf(beamNb).ray(rayNb).focusIx(bixelNb);

          if newIES % new IES
            voxNb = 1;
            iesNb = iesNb + 1;
            iesFocusIx = rayIESfocusIx;
            iesFocus = machine.data(energyIx).initFocus.SisFWHMAtIso(rayIESfocusIx);

            XMLplan(beamNb).IES(iesNb).energyIx = energyIx;
            XMLplan(beamNb).IES(iesNb).energy = iesEnergy;
            XMLplan(beamNb).IES(iesNb).focusIx = iesFocusIx;
            XMLplan(beamNb).IES(iesNb).focus = iesFocus;
            XMLplan(beamNb).IES(iesNb).bixels = [];

            newIES=false;

          else % check whether the focus from the new bixel is the same as the current focus
            iesFocusIxNewBixel = rayIESfocusIx;
            if ( iesFocusIx ~= iesFocusIxNewBixel)
              warndlg('ATTENTION: bixels in same IES with different foci!');
              warndlg('ATTENTION: only one focus is kept for consistence with the Syngo TPS at HIT!');
            end
          end % new IES

          XMLplan(beamNb).IES(iesNb).voxel(voxNb).x = voxel_x;
          XMLplan(beamNb).IES(iesNb).voxel(voxNb).y = voxel_y;
          XMLplan(beamNb).IES(iesNb).voxel(voxNb).particles = voxel_nbParticles;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,1) = voxel_x;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,2) = voxel_y;
          XMLplan(beamNb).IES(iesNb).voxelMatrix(voxNb,3) = voxel_nbParticles;
          XMLplan(beamNb).IES(iesNb).bixels(voxNb) = bixelIndex;

          XMLplan(beamNb).nbParticles = XMLplan(beamNb).nbParticles + voxel_nbParticles;
          voxNb = voxNb + 1;
        end % there are particles for given voxel
      elseif length(bixelNb)>1
        error('Unexpected number of IES in the same ray.');
        return;
      end % one IES found
    end % loop over rays
  end % find focus and "Voxel (x,y,nbParticles)" for each IES

  for iesNb = 1:length(XMLplan(beamNb).IES)

    spotNb = 1;

    voxelMatrix = sortrows(XMLplan(beamNb).IES(iesNb).voxelMatrix,[1,2]);

      while length(voxelMatrix)
          XMLplan(beamNb).IES(iesNb).rasterScanPath.x(spotNb,1) = voxelMatrix(1,1);
          XMLplan(beamNb).IES(iesNb).rasterScanPath.y(spotNb,1) = voxelMatrix(1,2);
          XMLplan(beamNb).IES(iesNb).rasterScanPath.particles(spotNb,1) = voxelMatrix(1,3);

          spotNb = spotNb + 1;

          tmpY = voxelMatrix(1,2);
          voxelMatrix(1,:) = [];

          if length(voxelMatrix)
            if (voxelMatrix(1,2) < tmpY)
              voxelMatrix = sortrows(voxelMatrix,[1,-2]);
            else
              voxelMatrix = sortrows(voxelMatrix,[1,2]);
            end
          end
      end
  end
end % loop over beams
