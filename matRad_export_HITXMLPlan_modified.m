%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Description: Function to export a matRad-based treatment plan in the HIT-XML
%               format based on the R-package HITXML from Steffen Greilich
%               available at https://r-forge.r-project.org/projects/hitxml/
%
%  Author:      Lucas Norberto Burigo
%
%  Please, report any problems to l.burigo@dkfz.de
%
%  Note1: the parameter minNbParticlesSpot is used to exclude spots with
%         number of particles below the threshold which can be delivered at HIT
%  Note2: the assigment of parameters for TxTable and BAMS needs to be verified
%  Note3: TxRoom fixed to Room3Gantry
%  Note4: parameters for Beam, Patient and TxInitiation assigned to dumb values
%  Note5: for execution with Octave, check the Xerces Java Parser path below
%
%  ScanPath: 'stfMode'  all spots are written in the plan file according to the
%  sorting of the stf 
%            'backforth' for a constant y the scan direction is forth and
%             for the next y back
%            'HIT' wäre schön
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function matRad_export_HITXMLPlan_modified(planFilename, pln, stf, resultGUI, scanPath, vis)

if nargin < 6
    vis = 0;
end

disp('HITXML exporter: Exporting plan in the HITXML format')

if ~strcmp(pln.radiationMode,'protons') && ~strcmp(pln.radiationMode,'carbon')
  error('HITXML plan for this radiationMode not supported!');
end
minNbParticlesSpot = 500000/1e6;  %for protons
minNrParticlesIES = 25000000;    %for protons
if strcmp(pln.radiationMode,'carbon')
     minNbParticlesSpot = 15000/1e6;   
     minNrParticlesIES = 0;
end


if ~strcmp(pln.machine,'HIT') 
  error('HITXML exporter only supports machine HIT');
end

% Compute bixel index as it is done on matRad_calParticleDose
counter = 0;

for i = 1:pln.numOfBeams; % loop over all beams

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

% OCTAVE dependencies
if exist('OCTAVE_VERSION','builtin');

  pkg load io;

  % check Xerces Java Parser dependencies
  [status cmdout] = system('if test -e /usr/share/java/xerces-j2.jar; then exit 0; echo exit 2; fi');
  if status == 2
   error('Xerces Java Parser dependencies not found');
    return;
  end

  [status cmdout] = system('if test -e /usr/share/java/xml-commons-apis.jar; then exit 0; echo exit 2; fi');
  if status == 2
   error('Xerces Java Parser dependencies not found');
    return;
  end

  javaaddpath ('/usr/share/java/xerces-j2.jar');
  javaaddpath ('/usr/share/java/xml-commons-apis.jar');
end

for beamNb = 1:pln.numOfBeams % loop over beams
  filename = sprintf('PBP_%02d_%s.xml',beamNb-1,planFilename);

  % helper function for energy selection
  round2 = @(a,b)round(a*10^b)/10^b;

if exist('OCTAVE_VERSION','builtin');
  % OCTAVE
  docNode = javaObject ('org.apache.xerces.dom.DocumentImpl');
  docNode.appendChild (docNode.createElement ('PTTxPlanMd5'));
else
  % MATLAB
  docNode = com.mathworks.xml.XMLUtils.createDocument('PTTxPlanMd5');
end

  docNode.appendChild(docNode.createComment('TREATMENT PLAN CREATED WITH MATRAD'));

  PTTxPlanMd5 = docNode.getDocumentElement;
  %PTTxPlanMd5.setAttribute('md5','noMD5');
  PTTxPlanMd5.setAttribute('xmlns:xsi','http://www.w3.org/2001/XMLSchema-instance');
  PTTxPlanMd5.setAttribute('xsi:noNamespaceSchemaLocation','RTT-PT-Plan.xsd');
  PTTxPlanMd5.setAttribute('md5','noMD5');

  PTTxPlan = docNode.createElement('PTTxPlan');
  PTTxPlanMd5.appendChild(PTTxPlan);

  Beam = docNode.createElement('Beam');
  Beam.setAttribute('uid','bee035c5-03f6-4e9c-94b9-31f0fc484db1');
  PTTxPlan.appendChild(Beam);

  RstFormat = docNode.createElement('RstFormat');
  RstFormat.appendChild(docNode.createTextNode('PT_2004'));
  Beam.appendChild(RstFormat);

  Patient = docNode.createElement('Patient');
  Patient.setAttribute('id','PT-2004-01');
  Patient.setAttribute('name','unknown');
  Patient.setAttribute('sex','unknown');
  Patient.setAttribute('birthDate','2001-01-01');
  Beam.appendChild(Patient);

  TxInitiation = docNode.createElement('TxInitiation');
  TxInitiation.setAttribute('therapist','None');
  TxInitiation.setAttribute('dateTime','2007-01-23T13:52:27.2343750+01:00');
  Beam.appendChild(TxInitiation);

  TxRoom = docNode.createElement('TxRoom');
%  TxRoom.setAttribute('name','Room1Fixed90');
  TxRoom.setAttribute('name','Room3');
  if strcmp(pln.radiationMode,'protons')
    TxRoom.setAttribute('projectile','PROTON');
    TxRoom.setAttribute('charge','');
    TxRoom.setAttribute('mass','');
    TxRoom.setAttribute('atomicNumber','');
  elseif strcmp(pln.radiationMode,'carbon')
    TxRoom.setAttribute('projectile','ION');
    TxRoom.setAttribute('charge','6');
    TxRoom.setAttribute('mass','12');
    TxRoom.setAttribute('atomicNumber','6');
  end
  Beam.appendChild(TxRoom);

  BAMS = docNode.createElement('BAMS');
  if strcmp(pln.radiationMode,'protons')
    BAMS.setAttribute('rippleFilter','254'); % value obtained from a xml plan sample from HIT
    BAMS.setAttribute('rangeShifter','254'); % value obtained from a xml plan sample from HIT
  elseif strcmp(pln.radiationMode,'carbon')
    BAMS.setAttribute('rippleFilter','3');   % value obtained from a xml plan sample from HIT
    BAMS.setAttribute('rangeShifter','254'); % value obtained from a xml plan sample from HIT
  end
  BAMS.setAttribute('rangeShifterDistance','0'); % value obtained from a xml plan sample from HIT
  Beam.appendChild(BAMS);

  TxTable = docNode.createElement('TxTable');
  TxTable.setAttribute('roll','0');         % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('pitch','0');        % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('lateral','0');      % value obtained from a xml plan sample from HIT
  TxTable.setAttribute('longitudinal','0'); % value obtained from a xml plan sample from HIT
  aValue = sprintf('%d',stf(beamNb).couchAngle);
  TxTable.setAttribute('isocentricAngle',aValue);
  TxTable.setAttribute('vertical','0');     % value obtained from a xml plan sample from HIT
  Beam.appendChild(TxTable);

  Gantry = docNode.createElement('Gantry');
  aValue = sprintf('%d',stf(beamNb).gantryAngle);
  Gantry.setAttribute('angle',aValue);
  Beam.appendChild(Gantry);

  % Find IES values
  iesArray = [];
  for rayNb=1:stf(beamNb).numOfRays
    iesArray = unique([iesArray stf(beamNb).ray(rayNb).energy]);
  end

  iesNb=0; % used to get rid off IES in the iesArray with no particles (weigth = 0)
  for iesArrayIx=1:length(iesArray) % find focus and "Voxel (x,y,nbParticles)" for each IES
    clear iesFocus;
    clear iesEnergy;
    newIES=false;
    iesEnergy = iesArray(iesArrayIx);
    iesFocus = 0;
    counter = 1;
    tmpIES = [];
    
    for rayNb=1:stf(beamNb).numOfRays % loop over rays

      % find index of used energy (round to keV for numerical reasons
      bixelNb = find(round2(stf(beamNb).ray(rayNb).energy,4) == round2(iesEnergy,4)==1);
      energyIx = find(round2(iesEnergy,4) == round2(availableEnergies,4)==1);

      if length(bixelNb)==1 % one IES found

        bixelIndex = find([tmp.beamNum==beamNb & tmp.rayNum==rayNb & tmp.bixelNum==bixelNb]==1);

        voxel_nbParticles = resultGUI.w(bixelIndex); % if minNrParticlesSpot is same for postprocessing and export_Plan resultGUI.w would give same result
        voxel_nbParticles = round(1e6*voxel_nbParticles);

        % check whether there are (enough) particles for beam delivery
        if (voxel_nbParticles>minNbParticlesSpot)
            
            rayPos_bev = stf(beamNb).ray(rayNb).rayPos_bev;
            % HITXML, x-y-z(beam)
            % matRad, x-y(beam)-z -> z-x-y(beam)
            voxel_x = rayPos_bev(3);
            voxel_y = rayPos_bev(1);
  
            rayIESenergy=stf(beamNb).ray(rayNb).energy(bixelNb);
            rayIESfocusIx=stf(beamNb).ray(rayNb).focusIx(bixelNb);
         
            % check whether the focus from the new bixel is the same as the current focus
            iesFocusIxNewBixel = rayIESfocusIx;
            if(~newIES)
                iesFocusIx = rayIESfocusIx;
            end
            if ( iesFocusIx ~= iesFocusIxNewBixel)
              warndlg('ATTENTION: bixels in same IES with different foci!');
              warndlg('ATTENTION: only one focus is kept for consistence with the Syngo TPS at HIT!');
            end
                
           tmpIES.voxel_x(counter) = voxel_x;
           tmpIES.voxel_y(counter) = voxel_y;
           tmpIES.voxel_nbParticles(counter) = voxel_nbParticles;
           counter= counter +1;
           
           newIES = true;
        end % there are enough particles            
      elseif length(bixelNb)>1
          error('Unexpected number of IES in the same ray.');
          return;
      end % one IES found
    end % loop over rays
    if newIES % new IES
        %check if enough particles in IES
        if(sum(tmpIES.voxel_nbParticles(:)) < minNrParticlesIES)
            fprintf(['Not enough particles in IES' num2str(rayIESenergy) 'deleted']);
        else    
        
        iesFocus = machine.data(energyIx).initFocus.SisFWHMAtIso(rayIESfocusIx);

        IES = docNode.createElement('IES');
        aValue = sprintf('%d',iesNb);
        IES.setAttribute('number',aValue); 
        aValue = sprintf('%.2f',iesEnergy);
        IES.setAttribute('energy',aValue);
        aValue = sprintf('%.1f',iesFocus);
        IES.setAttribute('focus',aValue);

         Beam.appendChild(IES);
         iesNb = iesNb+1;  
         
         if(strcmp(scanPath, 'stfMode'))
             for c= 1:counter-1
                 Voxel = docNode.createElement('Voxel');
                 aValue = sprintf('%f',tmpIES.voxel_x(c));
                 Voxel.setAttribute('x',aValue);
                 aValue = sprintf('%f',tmpIES.voxel_y(c));
                 Voxel.setAttribute('y',aValue);
                 aValue = sprintf('%d',tmpIES.voxel_nbParticles(c));
                 Voxel.setAttribute('particles',aValue);
                 IES.appendChild(Voxel);
             end
             
             if(vis)
                 figure
              numOfCities = counter-1;      
             cities = [tmpIES.voxel_x(:), tmpIES.voxel_y(:)].';
             hold on
             plot(cities(1,:),cities(2,:),'r.','MarkerSize',30)
             plot(cities(1,:),cities(2,:),'g','LineWidth',2)
                          title(['Energy slice ' num2str(iesEnergy) ' MeV']);    
             xlabel('x [mm]')
             ylabel('y [mm]')
             end
         elseif(strcmp(scanPath,'backforth'))
             c=1;
             voxel_y = tmpIES.voxel_y(c);
             while(c<counter)  
                 while(c<counter && tmpIES.voxel_y(c) == voxel_y)
                     Voxel = docNode.createElement('Voxel');
                     aValue = sprintf('%f',tmpIES.voxel_x(c));
                     Voxel.setAttribute('x',aValue);
                     aValue = sprintf('%f',tmpIES.voxel_y(c));
                     Voxel.setAttribute('y',aValue);
                     aValue = sprintf('%d',tmpIES.voxel_nbParticles(c));
                     Voxel.setAttribute('particles',aValue);
                     IES.appendChild(Voxel);
                     c=c+1;
                 end
                 if( c<counter)
                     voxel_y = tmpIES.voxel_y(c);
                 end
                 c_start = c;
                 while(c<counter && tmpIES.voxel_y(c) == voxel_y)
                     c=c+1;
                 end
                 c_stop = c-1;
                 for c= c_stop:-1:c_start
                     Voxel = docNode.createElement('Voxel');
                     aValue = sprintf('%f',tmpIES.voxel_x(c));
                     Voxel.setAttribute('x',aValue);
                     aValue = sprintf('%f',tmpIES.voxel_y(c));
                     Voxel.setAttribute('y',aValue);
                     aValue = sprintf('%d',tmpIES.voxel_nbParticles(c));
                     Voxel.setAttribute('particles',aValue);
                     IES.appendChild(Voxel);
                 end
                 c= c_stop +1;
                 if( c<counter)
                     voxel_y = tmpIES.voxel_y(c);
                 end
             end
             if(vis)
             cities = [];                
             c=1;
             n=1;
             voxel_y = tmpIES.voxel_y(c);
             while(c<counter)  
                 while(c<counter && tmpIES.voxel_y(c) == voxel_y)
                     cities(1,n) = tmpIES.voxel_x(c);
                     cities(2,n) = tmpIES.voxel_y(c);
                     n = n+1;
                     c= c+1;
                 end
             
             if( c<counter)
                     voxel_y = tmpIES.voxel_y(c);
                 end
                 c_start = c;
                 while(c<counter && tmpIES.voxel_y(c) == voxel_y)
                     c=c+1;
                 end
                 c_stop = c-1;
                 for c= c_stop:-1:c_start
                     cities(1,n) = tmpIES.voxel_x(c);
                     cities(2,n) = tmpIES.voxel_y(c);
                     n = n+1;
                 end
                 c= c_stop +1;
                 if( c<counter)
                     voxel_y = tmpIES.voxel_y(c);
                 end
             end
              figure
             hold on
             plot(cities(1,:),cities(2,:),'r.','MarkerSize',30)
             plot(cities(1,:),cities(2,:),'g','LineWidth',2)
             title(['Energy slice ' num2str(iesEnergy) ' MeV']);    
             xlabel('x [mm]')
             ylabel('y [mm]')

             
             end
         elseif(strcmp(scanPath,'TSP'))
             numOfCities = counter-1;      
             cities = [tmpIES.voxel_x(:), tmpIES.voxel_y(:)].';
             
             %genetic algorith for traveling salesman problem: 
             %cities = x,y positionen der spots in IES
             distances = zeros(numOfCities,numOfCities);
             for i = 1:numOfCities
                 for j = i+1:numOfCities
                     distances(i,j) = sum((cities(:,i)-cities(:,j)).^2)^.5;
                     distances(j,i) = distances(i,j);
                 end
             end
             clf
             hold on
             plot(cities(1,:),cities(2,:),'r.','MarkerSize',30)
             trajectoryVisHandle = [];
             
             numOfVisPoints = 2;
             visPoints = NaN*ones(numOfVisPoints,2,numOfCities,numOfCities);
             for i = 1:numOfCities
                 for j = i+1:numOfCities
                     visPoints(1,:,i,j) = cities(:,i);
                     visPoints(numOfVisPoints,:,i,j) = cities(:,j);
                     visPoints(:,:,j,i) = visPoints(:,:,i,j);
                 end
             end
             trajectory = NaN*ones(numOfVisPoints*(numOfCities-1),2);
    
             % set up evolutionary process
             numOfGenerations = 10000;
             objFunc          = NaN*ones(numOfGenerations,1);
             lag              = 100; % for convergence measure
             breakBool        = 0;
             numOfChromosomes = 1000;
             slideProb        = 1;
             swapProb         = 1;
             flipProb         = 1;
             cumSlideProb     = slideProb/sum([slideProb swapProb flipProb]);
             cumSwapProb      = (slideProb+swapProb)/sum([slideProb swapProb flipProb]);
             survivalRate     = .05;
             numOfSurvivors   = ceil(survivalRate*numOfChromosomes);

             population    = NaN*ones(numOfChromosomes,numOfCities);
             for i = 1:numOfChromosomes
                 population(i,:) = randperm(numOfCities);
             end

             generation = 0;
             while (~breakBool) && generation < numOfGenerations
                 generation = generation + 1;
                 
                 % assess fitness of all chromosomes
                 fitness = sum(distances((population(:,1:end-1)-1)*numOfCities + population(:,2:end)),2);

                 objFunc(generation) = min(fitness);
                 if generation > lag
                     breakBool = numel(unique(objFunc(generation-100:generation))) == 1; 
                 end
    
                 % sort population and fitness
                 [tmpTSP,ranking]   = sort(fitness);
                 population = population(ranking,:);

                 % genetic operators
                 randTriggers = rand(numOfChromosomes-numOfSurvivors,1);
                 randLoci     = ceil(numOfCities*rand(numOfChromosomes-numOfSurvivors,2));
                 for i = numOfSurvivors+1:numOfChromosomes
                     ix     = mod(i,numOfSurvivors)+1;
                     parent = population(ix,:);
                     prob   = randTriggers(i-numOfSurvivors);
                     if prob < cumSlideProb % slide
                         population(i,:) = parent([randLoci(i-numOfSurvivors,1):end 1:randLoci(i-numOfSurvivors,1)-1]);
                     elseif prob < cumSwapProb % swap
                         population(i,:) = parent;
                         population(i,randLoci(i-numOfSurvivors,1)) = parent(randLoci(i-numOfSurvivors,2));
                         population(i,randLoci(i-numOfSurvivors,2)) = parent(randLoci(i-numOfSurvivors,1));
                     else % flip
                         population(i,:) = parent;
                         if randLoci(i-numOfSurvivors,1) < randLoci(i-numOfSurvivors,2)
                         population(i,randLoci(i-numOfSurvivors,1):randLoci(i-numOfSurvivors,2)) = population(i,randLoci(i-numOfSurvivors,2):-1:randLoci(i-numOfSurvivors,1));
                         else
                         population(i,randLoci(i-numOfSurvivors,2):randLoci(i-numOfSurvivors,1)) = population(i,randLoci(i-numOfSurvivors,1):-1:randLoci(i-numOfSurvivors,2));                
                         end
                    end
                 end
    
                 % plot fittest individual    
                 delete(trajectoryVisHandle)
                 trajectoryVisHandle = ones(numOfCities-1,1);
                 for i = 1:numOfCities-1
                     start = population(1,i);
                     stop  = population(1,i+1);

                     trajectoryVisHandle(i) = plot(visPoints(:,1,start,stop),visPoints(:,2,start,stop),'g','LineWidth',2);
        
                 end
                 title(['Generation # ' num2str(generation)])
                 drawnow;
             end

             trajectory = cities(:,population(1,:));
        
             %write best trajectory in XML Plan file
             c = 1;
             while(c<counter)
                 Voxel = docNode.createElement('Voxel');
                 aValue = sprintf('%f',trajectory(1,c));
                 Voxel.setAttribute('x',aValue);
                 aValue = sprintf('%f',trajectory(2,c));
                 Voxel.setAttribute('y',aValue);
                 h = find(tmpIES.voxel_x == trajectory(1,c) & tmpIES.voxel_y == trajectory(2,c));
                 aValue = sprintf('%d',tmpIES.voxel_nbParticles(h));
                 Voxel.setAttribute('particles',aValue);
                 IES.appendChild(Voxel);
                 c=c+1;
             end
             
         else                             
             error('No available Scan Path Mode selected');             
         end
        end
    end   
    
  end % find focus and "Voxel (x,y,nbParticles)" for each IES
  
  xmlwrite(filename,docNode);
end % loop over beams

disp('HITXML exporter: End of exporter function.')
