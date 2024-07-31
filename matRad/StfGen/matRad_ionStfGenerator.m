classdef matRad_ionStfGenerator < matRad_externalStfGenerator     

        properties (Constant)
        name = 'ionStfGen';
        shortName = 'ionStfGen';
    end 
    
    
    properties (Access = protected)
        availableEnergies
        pbMargin
        availablePeakPosRaShi
        longitudinalSpotSpacing

    end



    methods 
        function this = matRad_ionStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_externalStfGenerator(pln);
            this.radiationMode = {'protons','carbon','helium'};

            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));
  



            this.availableEnergies = availableEnergies;
            this.pbMargin = pbMargin;
            this.availablePeakPosRaShi = availablePeakPosRaShi;
            this.longitudinalSpotSpacing = longitudinalSpotSpacing;

            
            if ~isfield(pln, 'propStf')
                this.matRad_cfg.dispError('no applicator information in pln struct');
            end
         end
    end

    methods 

        function stf = generate(this, ct, cst)
            stf = generate@matRad_externalStfGenerator(this, ct, cst);
            this.initializePatientGeometry(ct, cst);
            stf = this.generateSourceGeometry(ct, cst);
        end
    end

    methods (Access = protected)
        function initializePatientGeometry(this, ct, cst)
            % Initialize the patient geometry
            initializePatientGeometry@matRad_externalStfGenerator(this);

            if ~isfield(pln.propStf,'useRangeShifter')
            pln.propStf.useRangeShifter = false;
            end

            this.availableEnergies = [machine.data.energy];
            availablePeakPos  = [machine.data.peakPos] + [machine.data.offset];
            availableWidths   = [machine.data.initFocus];
            availableWidths   = [availableWidths.SisFWHMAtIso];
            maxPBwidth        = max(availableWidths) / 2.355;

            if pln.propStf.useRangeShifter
            %For now only a generic range shifter is used whose thickness is
            %determined by the minimum peak width to play with
            rangeShifterEqD = round(min(availablePeakPos)* 1.25);
            this.availablePeakPosRaShi = availablePeakPos - rangeShifterEqD;

            matRad_cfg.dispWarning('Use of range shifter enabled. matRad will generate a generic range shifter with WEPL %f to enable ranges below the shortest base data entry.',rangeShifterEqD);
            end

            if ~isfield(pln.propStf, 'longitudinalSpotSpacing')
                this.longitudinalSpotSpacing = matRad_cfg.propStf.defaultLongitudinalSpotSpacing;
            else
                this.longitudinalSpotSpacing = pln.propStf.longitudinalSpotSpacing;
            end

            if sum(availablePeakPos<0)>0
                matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file')
            end

        end
        
        function pbMargin = getPbMargin(this)
            pbMargin = getPbMargin@matRad_externalStfGenerator(this);
            %Compute a margin to account for pencil beam width
            pbMargin = min(maxPBwidth,pln.propStf.bixelWidth);
        end
        
        function stf = generateSourceGeometry(this, ct, cst)
            stf = generateSourceGeometry@matRad_externalStfGenerator(this, ct, cst);

        end

        function initializeEnergy(this,ct,cst)
            initializeEnergy@matRad_externalStfGenerator(this);

                % target hit
                rhoVOITarget = [];
                for shiftScen = 1:pln.multScen.totNumShiftScen
                    rhoVOITarget = [rhoVOITarget, rho{shiftScen}{end}];
                end

                if any(rhoVOITarget)
                    Counter = 0;

                    %Here we iterate through scenarios to check the required
                    %energies w.r.t lateral position.
                    %TODO: iterate over the linear scenario mask instead?
                    for CtScen = 1:pln.multScen.numOfCtScen
                        for ShiftScen = 1:pln.multScen.totNumShiftScen
                            for RangeShiftScen = 1:pln.multScen.totNumRangeScen
                                if pln.multScen.scenMask(CtScen,ShiftScen,RangeShiftScen)
                                    Counter = Counter+1;

                                    % compute radiological depths
                                    % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                                    radDepths = cumsum(l{ShiftScen} .* rho{ShiftScen}{CtScen});

                                    if pln.multScen.relRangeShift(RangeShiftScen) ~= 0 || pln.multScen.absRangeShift(RangeShiftScen) ~= 0
                                        radDepths = radDepths +...                                                        % original cube
                                            rho{ShiftScen}{CtScen}*pln.multScen.relRangeShift(RangeShiftScen) +... % rel range shift
                                            pln.multScen.absRangeShift(RangeShiftScen);                           % absolute range shift
                                        radDepths(radDepths < 0) = 0;
                                    end

                                    % find target entry & exit
                                    diff_voi    = [diff([rho{ShiftScen}{end}])];
                                    entryIx = find(diff_voi == 1);
                                    exitIx = find(diff_voi == -1);

                                    %We approximate the interface using the
                                    %rad depth between the last voxel before
                                    %and the first voxel after the interface
                                    %This captures the case that the first
                                    %relevant voxel is a target voxel
                                    targetEntry(Counter,1:length(entryIx)) = (radDepths(entryIx) + radDepths(entryIx+1)) ./ 2;
                                    targetExit(Counter,1:length(exitIx)) = (radDepths(exitIx) + radDepths(exitIx+1)) ./ 2;
                                end
                            end
                        end
                    end

                    targetEntry(targetEntry == 0) = NaN;
                    targetExit(targetExit == 0)   = NaN;

                    targetEntry = min(targetEntry);
                    targetExit  = max(targetExit);

                    %check that each energy appears only once in stf
                    if(numel(targetEntry)>1)
                        m = numel(targetEntry);
                        while(m>1)
                            if(targetEntry(m) < targetExit(m-1))
                                targetExit(m-1) = max(targetExit(m-1:m));
                                targetExit(m)=[];
                                targetEntry(m-1) = min(targetEntry(m-1:m));
                                targetEntry(m)=[];
                                m = numel(targetEntry)+1;
                            end
                            m=m-1;
                        end
                    end

                    if numel(targetEntry) ~= numel(targetExit)
                        matRad_cfg.dispError('Inconsistency during ray tracing. Please check correct assignment and overlap priorities of structure types OAR & TARGET.');
                    end

                    stf(i).ray(j).energy = [];
                    stf(i).ray(j).rangeShifter = [];

                    % Save energies in stf struct
                    for k = 1:numel(targetEntry)

                        %If we need lower energies than available, consider
                        %range shifter (if asked for)
                        if any(targetEntry < min(availablePeakPos)) && pln.propStf.useRangeShifter
                            %Get Energies to use with range shifter to fill up
                            %non-reachable low-range spots
                            raShiEnergies = this.availableEnergies(this.availablePeakPosRaShi >= targetEntry(k) & min(availablePeakPos) > this.availablePeakPosRaShi);

                            raShi.ID = 1;
                            raShi.eqThickness = rangeShifterEqD;
                            raShi.sourceRashiDistance = round(ctEntryPoint - 2*rangeShifterEqD,-1); %place a little away from entry, round to cms to reduce number of unique settings

                            stf(i).ray(j).energy = [stf(i).ray(j).energy raShiEnergies];
                            stf(i).ray(j).rangeShifter = [stf(i).ray(j).rangeShifter repmat(raShi,1,length(raShiEnergies))];
                        end

                        %Normal placement without rangeshifter
                        newEnergies = this.availableEnergies(availablePeakPos>=targetEntry(k)&availablePeakPos<=targetExit(k));


                        stf(i).ray(j).energy = [stf(i).ray(j).energy newEnergies];


                        raShi.ID = 0;
                        raShi.eqThickness = 0;
                        raShi.sourceRashiDistance = 0;
                        stf(i).ray(j).rangeShifter = [stf(i).ray(j).rangeShifter repmat(raShi,1,length(newEnergies))];
                    end


                    targetEntry = [];
                    targetExit = [];


                    % book keeping & calculate focus index
                    stf(i).numOfBixelsPerRay(j) = numel([stf(i).ray(j).energy]);
                    currentMinimumFWHM = matRad_interp1(machine.meta.LUT_bxWidthminFWHM(1,:)',...
                        machine.meta.LUT_bxWidthminFWHM(2,:)',...
                        pln.propStf.bixelWidth, ...
                        machine.meta.LUT_bxWidthminFWHM(2,end));
                    focusIx  =  ones(stf(i).numOfBixelsPerRay(j),1);
                    [~, vEnergyIx] = min(abs(bsxfun(@minus,[machine.data.energy]',...
                        repmat(stf(i).ray(j).energy,length([machine.data]),1))));

                    % get for each spot the focus index
                    for k = 1:stf(i).numOfBixelsPerRay(j)
                        focusIx(k) = find(machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                    end

                    stf(i).ray(j).focusIx = focusIx';

                    %Get machine bounds
                    numParticlesPerMU = 1e6*ones(1,stf(i).numOfBixelsPerRay(j));
                    minMU = zeros(1,stf(i).numOfBixelsPerRay(j));
                    maxMU = Inf(1,stf(i).numOfBixelsPerRay(j));
                    for k = 1:stf(i).numOfBixelsPerRay(j)
                        if isfield(machine.data(vEnergyIx(k)),'MUdata')
                            MUdata = machine.data(vEnergyIx(k)).MUdata;
                            if isfield(MUdata,'numParticlesPerMU')
                                numParticlesPerMU(k) = MUdata.numParticlesPerMU;
                            end

                            if isfield(MUdata,'minMU')
                                minMU(k) = MUdata.minMU;
                            end

                            if isfield(MUdata,'maxMU')
                                maxMU(k) = MUdata.maxMU;
                            end
                        end
                    end

                    stf(i).ray(j).numParticlesPerMU = numParticlesPerMU;
                    stf(i).ray(j).minMU = minMU;
                    stf(i).ray(j).maxMU = maxMU;

                else % target not hit
                    stf(i).ray(j)               = [];
                    stf(i).numOfBixelsPerRay(j) = [];
                end


        end

        function initializeIonPostProcessing(this,ct,cst)
            initializeIonPostProcessing@matRad_externalStfGenerator(this);

            % get minimum energy per field
            minEnergy = min([stf(i).ray.energy]);
            maxEnergy = max([stf(i).ray.energy]);

            % get corresponding peak position
            minPeakPos  = machine.data(minEnergy == this.availableEnergies).peakPos;
            maxPeakPos  = machine.data(maxEnergy == this.availableEnergies).peakPos;

            % find set of energyies with adequate spacing


            stf(i).longitudinalSpotSpacing = this.longitudinalSpotSpacing;

            tolerance              = this.longitudinalSpotSpacing/10;

            useEnergyBool = availablePeakPos >= minPeakPos & availablePeakPos <= maxPeakPos;

            ixCurr = find(useEnergyBool,1,'first');
            ixRun  = ixCurr + 1;
            ixEnd  = find(useEnergyBool,1,'last');

            while ixRun <= ixEnd
                if abs(availablePeakPos(ixRun)-availablePeakPos(ixCurr)) < ...
                        this.longitudinalSpotSpacing - tolerance
                    useEnergyBool(ixRun) = 0;
                else
                    ixCurr = ixRun;
                end
                ixRun = ixRun + 1;
            end

            for j = stf(i).numOfRays:-1:1
                for k = stf(i).numOfBixelsPerRay(j):-1:1
                    maskEnergy = stf(i).ray(j).energy(k) == this.availableEnergies;
                    if ~useEnergyBool(maskEnergy)
                        stf(i).ray(j).energy(k)         = [];
                        stf(i).ray(j).focusIx(k)        = [];
                        stf(i).ray(j).rangeShifter(k)   = [];
                        stf(i).numOfBixelsPerRay(j) = stf(i).numOfBixelsPerRay(j) - 1;
                    end
                end
                if isempty(stf(i).ray(j).energy)
                    stf(i).ray(j) = [];
                    stf(i).numOfBixelsPerRay(j) = [];
                    stf(i).numOfRays = stf(i).numOfRays - 1;
                end
            end

        end
    end
end
