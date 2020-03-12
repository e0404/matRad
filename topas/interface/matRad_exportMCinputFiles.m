function matRad_exportMCinputFiles(pln, ct, stf, MCsettings, w)     % Script to export MC input files

if ~isfield(MCsettings,'MC_PBS')
    MCsettings.MC_PBS = 0;
end

if ~exist('w','var')
    % default: weights set to one
    w = ones(stf.totalNumOfBixels,1);
end

if ~isfield(MCsettings,'ElectronProdCut_mm')
    MCsettings.ElectronProdCut_mm = 0.5;
end

if MCsettings.MC_PBS
  disp('Using PBS mode');
else
  disp('Not using PBS mode. Far SAD mode!');
end

nbRuns = 5;

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

energySpectrum = [];
hasBeamEnergySpectrum = false;
minNbSpectrumPoints = Inf;
maxNbSpectrumPoints = 0;

% load machine file
fileName = [pln.radiationMode '_' pln.machine];
try
  load(fileName);
catch
  error(['Could not find the following machine file: ' fileName ]);
end

% NozzleAxialDistance
nozzleAxialDistance_mm = 1500;
SAD_mm  = machine.meta.SAD;
if isfield(machine.meta,'nozzleAxialDistance')
  disp('Using NAD from basedata')
  nozzleAxialDistance_mm = machine.meta.nozzleAxialDistance;
else
  disp('Using default nozzleAxialDistance')
end

if isfield(MCsettings,'MCcode') && strcmp('TOPAS',MCsettings.MCcode)
  disp('Checking basedata for beam energy spectrum...');
  if isfield(machine.data,'energySpectrum')
    hasBeamEnergySpectrum = true;
    disp('Beam energy spectrum available');
    for ies=1:length(machine.data)
      minNbSpectrumPoints = min([length(machine.data(ies).energySpectrum.energy_MeVpN),minNbSpectrumPoints]);
      maxNbSpectrumPoints = max([length(machine.data(ies).energySpectrum.energy_MeVpN),maxNbSpectrumPoints]);
      energySpectrum = [energySpectrum;machine.data(ies).energySpectrum];
    end
    if minNbSpectrumPoints == 1
      if maxNbSpectrumPoints > 1
        error('Data base has mixed monoenergetic/polienergetic IES. Aborting...')
      else
        energySpectrum = [];
        hasBeamEnergySpectrum = false;
        warning('Data base has IES for monoenergetic beams')    
      end
    end
  end
end

MCheader = [];
MCheader.nbRuns = nbRuns;

numOfBeams = [];
if isfield(pln,'propStf')
numOfBeams = pln.propStf.numOfBeams;
else
numOfBeams = pln.numOfBeams;
end
MCheader.nbFields = numOfBeams;
MCheader.tallies = {'physicalDose'};
MCheader.simLabel = 'matRad_plan';
MCheader.cubeDim = ct.cubeDim;
MCheader.RSP = ct.cube{1};
MCheader.voxelDimensions = ct.resolution;

XMLplan = matRad_create_XMLplan(pln,stf,w,MCsettings);

nbParticlesTotal = sum([XMLplan(:).nbParticles]);
%disp('Total number of particles for the plan: %s\n',num2str(nbParticlesTotal))

nbHistoriesTotal = uint32(MCsettings.fractionHistories*double(nbParticlesTotal));
%printf('Simulating %d histories in total\n',nbHistoriesTotal)

for beamNb = 1:numOfBeams % loop over beams

  nbParticles = XMLplan(beamNb).nbParticles;
  nbHistories = uint32(MCsettings.fractionHistories*double(nbParticles));
  %printf('Field %d with %s particles\n',beamNb,num2str(nbParticles))
  %printf('Simulating %s histories for field %d\n',num2str(nbHistories),beamNb)

  MCheader.nbHistories(beamNb) = nbHistories;
  MCheader.nbParticles(beamNb) = nbParticles;

  h = fopen([MCsettings.runsPath 'beamSetup_matRad_plan_field' num2str(beamNb) '.txt'],'w');

  fprintf(h,'i:Ts/ShowHistoryCountAtInterval = %d\n', nbHistories/20);
  fprintf(h,'s:Sim/PlanLabel = "simData_matRad_plan_field%d_run" + Ts/Seed\n', beamNb);
if isfield(pln,'propStf')
  fprintf(h,'d:Sim/GantryAngle = %f deg\n', pln.propStf.gantryAngles(beamNb));
  fprintf(h,'d:Sim/CouchAngle = %f deg\n', pln.propStf.couchAngles(beamNb));
else
  fprintf(h,'d:Sim/GantryAngle = %f deg\n', pln.gantryAngles(beamNb));
  fprintf(h,'d:Sim/CouchAngle = %f deg\n', pln.couchAngles(beamNb));
end
  % uses generic beamSetup

    if strcmp(XMLplan(beamNb).ion, '12C')
      fprintf(h,'s:Sim/ParticleName = "GenericIon(6,12,6)"\n');
      fprintf(h,'u:Sim/ParticleMass = 12.0\n');
    elseif strcmp(XMLplan(beamNb).ion, '1H')
      fprintf(h,'s:Sim/ParticleName = "proton"\n');
      fprintf(h,'u:Sim/ParticleMass = 1.0\n');
    end

    ies = [];
    energy = [];
    focus = [];
    posX = [];
    posY = [];
    current = [];
    for iesNb = 1:size(XMLplan(beamNb).IES,2)
      nbRasterPoints = size(XMLplan(beamNb).IES(iesNb).rasterScanPath.x,1);
      ies     = [ies;XMLplan(beamNb).IES(iesNb).energyIx*ones(nbRasterPoints,1)];
      energy  = [energy;XMLplan(beamNb).IES(iesNb).energy*ones(nbRasterPoints,1)];
      focus   = [focus;XMLplan(beamNb).IES(iesNb).focus*ones(nbRasterPoints,1)];
      posX    = [posX;-1.*XMLplan(beamNb).IES(iesNb).rasterScanPath.x];
      posY    = [posY;XMLplan(beamNb).IES(iesNb).rasterScanPath.y];
      current = [current;uint32(MCsettings.fractionHistories*double(XMLplan(beamNb).IES(iesNb).rasterScanPath.particles))];
    end

    idx = find(current < 1);
    energy(idx) = [];
    focus(idx) = [];
    posX(idx) = [];
    posY(idx) = [];
    current(idx) = [];

    % printf('Number of histories to simulate: %d\n',nbHistories);
    % printf('Total current: %d\n',sum(current));
    while sum(current) ~= nbHistories
      % Randomly pick an index with the weigth given by the current
      idx = 1:length(current);
      [~,R] = histc(rand(1),cumsum([0;double(current(:))./double(sum(current))]));
      randIx = idx(R);

      if (sum(current) > nbHistories)
        if current(randIx)>1
          current(randIx) = current(randIx)-1;
        end
      else
        current(randIx) = current(randIx)+1;
      end
      % printf('Adjusting total current to : %d\n',sum(current));
    end

    fprintf(h,'i:Sim/NbThreads = 0\n');

    fprintf(h,'d:Tf/TimelineStart = 0. ms\n');
    fprintf(h,'d:Tf/TimelineEnd = %g ms\n',10.*(length(posX)));
    fprintf(h,'i:Tf/NumberOfSequentialTimes = %d\n',length(posX));
    fprintf(h,'dv:Tf/Beam/Spot/Times = %d %s ms\n',length(posX),strtrim(sprintf('%.1f ',10.*(1:length(posX)))));

    fprintf(h,'s:Tf/Beam/Energy/Function = "Step"\n');
    fprintf(h,'dv:Tf/Beam/Energy/Times = Tf/Beam/Spot/Times ms\n');
    fprintf(h,'dv:Tf/Beam/Energy/Values = %d %s MeV\n',length(energy),strtrim(sprintf('%f ',energy)));

    fprintf(h,'s:Tf/Beam/FocusFWHM/Function = "Step"\n');
    fprintf(h,'dv:Tf/Beam/FocusFWHM/Times = Tf/Beam/Spot/Times ms\n');
    fprintf(h,'dv:Tf/Beam/FocusFWHM/Values = %d %s mm\n',length(focus),strtrim(sprintf('%f ',focus)));

    if MCsettings.MC_PBS
        % angleX corresponds to the rotation around the X axis necessary to move the spot in the Y direction
        % angleY corresponds to the rotation around the Y' axis necessary to move the spot in the X direction
        % note that Y' corresponds to the Y axis after the rotation of angleX around X axis
        angleX = atan(posY / SAD_mm);
        angleY = atan(-posX ./ (SAD_mm ./ cos(angleX)));
        posX = (posX / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);
        posY = (posY / SAD_mm)*(SAD_mm-nozzleAxialDistance_mm);

        fprintf(h,'s:Tf/Beam/AngleX/Function = "Step"\n');
        fprintf(h,'dv:Tf/Beam/AngleX/Times = Tf/Beam/Spot/Times ms\n');
        fprintf(h,'dv:Tf/Beam/AngleX/Values = %d %s rad\n',length(angleX),strtrim(sprintf('%f ',angleX)));

        fprintf(h,'s:Tf/Beam/AngleY/Function = "Step"\n');
        fprintf(h,'dv:Tf/Beam/AngleY/Times = Tf/Beam/Spot/Times ms\n');
        fprintf(h,'dv:Tf/Beam/AngleY/Values = %d %s rad\n',length(angleY),strtrim(sprintf('%f ',angleY)));
    end

    fprintf(h,'s:Tf/Beam/PosX/Function = "Step"\n');
    fprintf(h,'dv:Tf/Beam/PosX/Times = Tf/Beam/Spot/Times ms\n');
    fprintf(h,'dv:Tf/Beam/PosX/Values = %d %s mm\n',length(posX),strtrim(sprintf('%f ',posX)));

    fprintf(h,'s:Tf/Beam/PosY/Function = "Step"\n');
    fprintf(h,'dv:Tf/Beam/PosY/Times = Tf/Beam/Spot/Times ms\n');
    fprintf(h,'dv:Tf/Beam/PosY/Values = %d %s mm\n',length(posY),strtrim(sprintf('%f ',posY)));

    fprintf(h,'s:Tf/Beam/Current/Function = "Step"\n');
    fprintf(h,'dv:Tf/Beam/Current/Times = Tf/Beam/Spot/Times ms\n');
    fprintf(h,'iv:Tf/Beam/Current/Values = %d %s\n',length(current),strtrim(sprintf('%d ',current)));

    if hasBeamEnergySpectrum
      nbSpectrumPoints = length(machine.data(1).energySpectrum.energy_MeVpN);

      energyIx = zeros(1,length(energy));
      for ix=1:length(energy)
        energyIx(ix) = find(round2(energy(ix),4) == round2([machine.data.energy],4));
      end

      fprintf(h,'s:So/PencilBeam/BeamEnergySpectrumType = "Continuous"\n');
      fprintf(h,'dv:So/PencilBeam/BeamEnergySpectrumValues = %d %s MeV * Sim/ParticleMass\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Energy/Point%03d/Value ',1:nbSpectrumPoints)));
      fprintf(h,'uv:So/PencilBeam/BeamEnergySpectrumWeights = %d %s\n',nbSpectrumPoints,strtrim(sprintf('Tf/Beam/EnergySpectrum/Weight/Point%03d/Value ',1:nbSpectrumPoints)));

      points_energy = reshape([energySpectrum(energyIx).energy_MeVpN],[],length(energyIx));
      points_weight = reshape([energySpectrum(energyIx).weight],[],length(energyIx));
      for spectrumPoint=1:nbSpectrumPoints
        fprintf(h,'s:Tf/Beam/EnergySpectrum/Energy/Point%03d/Function = "Step"\n',spectrumPoint);
        fprintf(h,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
        fprintf(h,'dv:Tf/Beam/EnergySpectrum/Energy/Point%03d/Values = %d %s MeV\n',spectrumPoint,length(energy),strtrim(sprintf('%f ',points_energy(spectrumPoint,:))));
        fprintf(h,'s:Tf/Beam/EnergySpectrum/Weight/Point%03d/Function = "Step"\n',spectrumPoint);
        fprintf(h,'dv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Times = Tf/Beam/Spot/Times ms\n',spectrumPoint);
        fprintf(h,'uv:Tf/Beam/EnergySpectrum/Weight/Point%03d/Values = %d %s\n',spectrumPoint,length(energy),strtrim(sprintf('%f ',points_weight(spectrumPoint,:))));
      end
    else
      fprintf(h,'d:So/PencilBeam/BeamEnergy = Tf/Beam/Energy/Value MeV * Sim/ParticleMass\n');
    end

  MCheader.tallies = {'physicalDose'};

if isfield(pln,'propStf')
  fprintf(h,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-pln.propStf.isoCenter(beamNb,1));
  fprintf(h,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-pln.propStf.isoCenter(beamNb,2));
  fprintf(h,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-pln.propStf.isoCenter(beamNb,3));
else
  fprintf(h,'d:Ge/Patient/TransX      = %f mm\n',0.5*ct.resolution.x*(ct.cubeDim(2)+1)-pln.isoCenter(beamNb,1));
  fprintf(h,'d:Ge/Patient/TransY      = %f mm\n',0.5*ct.resolution.y*(ct.cubeDim(1)+1)-pln.isoCenter(beamNb,2));
  fprintf(h,'d:Ge/Patient/TransZ      = %f mm\n',0.5*ct.resolution.z*(ct.cubeDim(3)+1)-pln.isoCenter(beamNb,3));
end
  fprintf(h,'d:Ge/Patient/RotX=0. deg\n');
  fprintf(h,'d:Ge/Patient/RotY=0. deg\n');
  fprintf(h,'d:Ge/Patient/RotZ=0. deg\n');

  fprintf(h,'includeFile = %s\n', fullfile('.\matRad_RSPcube.txt'));
  fprintf(h,'###################\n\n');

  fbase = 0;

  % NozzleAxialDistance
  fprintf(h,'d:Ge/Nozzle/TransZ = -%f mm\n', nozzleAxialDistance_mm);
  if MCsettings.MC_PBS
      fprintf(h,'d:Ge/Nozzle/RotX = Tf/Beam/AngleX/Value rad\n');
      fprintf(h,'d:Ge/Nozzle/RotY = Tf/Beam/AngleY/Value rad\n');
      fprintf(h,'d:Ge/Nozzle/RotZ = 0.0 rad\n');
  end

  fprintf(h,['d:Ph/Default/CutForElectron = ',num2str(MCsettings.ElectronProdCut_mm),' mm\n']);

  fbase = fopen(['forwardMC_TOPAS_setup_generic_' pln.radiationMode '.txt'],'r');

  while ~feof(fbase)
    strLine = fgets(fbase); %# read line by line
    fprintf(h,'%s',strLine);
  end
  fclose(fbase);

  fclose(h);

  for runNb = 1:nbRuns
    h = fopen([MCsettings.runsPath 'forwardMC_matRad_plan_field' num2str(beamNb) '_run' num2str(runNb) '.txt'],'w');
    fprintf(h,'i:Ts/Seed = %d\n', runNb);
    fprintf(h,'includeFile = %s\n', fullfile(MCsettings.runsPath,['beamSetup_matRad_plan_field' num2str(beamNb) '.txt']));
    fclose(h);
  end % loop over runs
end % loop over beams

if exist('OCTAVE_VERSION','builtin');
  % OCTAVE
  save('-v7','MCparam.mat','MCheader');
else
  % MATLAB
  save('MCparam.mat','MCheader');
end

basematerial = '';
if ~exist('machine') || ~isfield(machine.meta,'basematerial')
  warning('Base material not defined in base data. Using Water')
  basematerial = 'Water';
else
  basematerial = machine.meta.basematerial;
end

matRad_exportCtTOPAS(ct, MCsettings.runsPath, basematerial)
end
