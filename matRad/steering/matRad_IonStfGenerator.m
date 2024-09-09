classdef matRad_IonStfGenerator < matRad_ExternalStfGenerator     

    properties (Constant)
        name = 'ionStfGen';
        shortName = 'ionStfGen';
        possibleRadiationModes = {'protons','helium','carbon'};
    end 

    properties
        longitudinalSpotSpacing;
        useRangeShifter = false;
    end
    
    properties (Access = protected)
        availableEnergies
        availablePeakPos
        availablePeakPosRaShi
        maxPBwidth
        pbMargin
    end



    methods 
        function this = matRad_IonStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_ExternalStfGenerator(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'protons';
            end
         end
    end


    methods (Access = protected)
        function initializePatientGeometry(this,ct, cst)
            % Initialize the patient geometry
            initializePatientGeometry@matRad_ExternalStfGenerator(this,ct, cst)
            matRad_cfg = MatRad_Config.instance;

            this.availableEnergies = [this.machine.data.energy];
            this.availablePeakPos  = [this.machine.data.peakPos] + [this.machine.data.offset];
            availableWidths   = [this.machine.data.initFocus];
            availableWidths   = [availableWidths.SisFWHMAtIso];
            this.maxPBwidth        = max(availableWidths) / 2.355;

            if this.useRangeShifter
                %For now only a generic range shifter is used whose thickness is
                %determined by the minimum peak width to play with
                rangeShifterEqD = round(min(this.availablePeakPos)* 1.25);
                this.availablePeakPosRaShi = this.availablePeakPos - rangeShifterEqD;

                matRad_cfg.dispWarning('Use of range shifter enabled. matRad will generate a generic range shifter with WEPL %f to enable ranges below the shortest base data entry.',rangeShifterEqD);
            end

            if sum(this.availablePeakPos<0)>0
                matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file')
            end

        end
        
        function pbMargin = getPbMargin(this)
            %Compute a margin to account for pencil beam width
            pbMargin = min(this.maxPBwidth,this.bixelWidth);
        end
        
        function stfElement = initializeEnergy(this,stfElement,ct) 

            isoCenterInCubeCoords = matRad_world2cubeCoords(stfElement.isoCenter,ct);

            for j = stfElement.numOfRays:-1:1

                for shiftScen = 1:this.multScen.totNumShiftScen
                        % ray tracing necessary to determine depth of the target
                        [alphas,l{shiftScen},rho{shiftScen},d12,~] = matRad_siddonRayTracer(isoCenterInCubeCoords + this.multScen.isoShift(shiftScen,:), ...
                            ct.resolution, ...
                            stfElement.sourcePoint, ...
                            stfElement.ray(j).targetPoint, ...
                            [ct.cube {this.voiTarget}]);

                        %Used for generic range-shifter placement
                        this.ctEntryPoint = alphas(1) * d12;
                end

                % target hit
                rhoVOITarget = [];
                for shiftScen = 1:this.multScen.totNumShiftScen
                    rhoVOITarget = [rhoVOITarget, rho{shiftScen}{end}];
                end

                if any(rhoVOITarget)
                    counter = 0;

                    %Here we iterate through scenarios to check the required
                    %energies w.r.t lateral position.
                    %TODO: iterate over the linear scenario mask instead?
                    for ctScen = 1:this.multScen.numOfCtScen
                        for shiftScen = 1:this.multScen.totNumShiftScen
                            for rangeShiftScen = 1:this.multScen.totNumRangeScen
                                if this.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                                    counter = counter+1;

                                    % compute radiological depths
                                    % http://www.ncbi.nlm.nih.gov/pubmed/4000088, eq 14
                                    radDepths = cumsum(l{shiftScen} .* rho{shiftScen}{ctScen});

                                    if this.multScen.relRangeShift(rangeShiftScen) ~= 0 || this.multScen.absRangeShift(rangeShiftScen) ~= 0
                                        radDepths = radDepths +...                                                      % original cube
                                            rho{shiftScen}{ctScen}*this.multScen.relRangeShift(rangeShiftScen) +...     % rel range shift
                                            this.multScen.absRangeShift(rangeShiftScen);                                % absolute range shift
                                        radDepths(radDepths < 0) = 0;
                                    end

                                    % find target entry & exit
                                    diff_voi    = [diff([rho{shiftScen}{end}])];
                                    entryIx = find(diff_voi == 1);
                                    exitIx = find(diff_voi == -1);

                                    %We approximate the interface using the rad depth between the last voxel before and the first voxel after the interface 
                                    % This captures the case that the first relevant voxel is a target voxel
                                    targetEntry(counter,1:length(entryIx))  = (radDepths(entryIx) + radDepths(entryIx+1)) ./ 2;
                                    targetExit(counter,1:length(exitIx))    = (radDepths(exitIx) + radDepths(exitIx+1))   ./ 2;
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

                    stfElement.ray(j).energy = [];
                    stfElement.ray(j).rangeShifter = [];

                    % Save energies in stf struct
                    for k = 1:numel(targetEntry)

                        %If we need lower energies than available, consider
                        %range shifter (if asked for)
                        if any(targetEntry < min(this.availablePeakPos)) && this.useRangeShifter
                            %Get Energies to use with range shifter to fill up
                            %non-reachable low-range spots
                            raShiEnergies = this.availableEnergies(this.availablePeakPosRaShi >= targetEntry(k) & min(this.availablePeakPos) > this.availablePeakPosRaShi);

                            raShi.ID = 1;
                            raShi.eqThickness = rangeShifterEqD;
                            raShi.sourceRashiDistance = round(ctEntryPoint - 2*rangeShifterEqD,-1); %place a little away from entry, round to cms to reduce number of unique settings

                            stfElement.ray(j).energy = [stfElement.ray(j).energy raShiEnergies];
                            stfElement.ray(j).rangeShifter = [stfElement.ray(j).rangeShifter repmat(raShi,1,length(raShiEnergies))];
                        end

                        %Normal placement without rangeshifter
                        newEnergies = this.availableEnergies(this.availablePeakPos>=targetEntry(k)&this.availablePeakPos<=targetExit(k));


                        stfElement.ray(j).energy = [stfElement.ray(j).energy newEnergies];


                        raShi.ID = 0;
                        raShi.eqThickness = 0;
                        raShi.sourceRashiDistance = 0;
                        stfElement.ray(j).rangeShifter = [stfElement.ray(j).rangeShifter repmat(raShi,1,length(newEnergies))];
                    end


                    targetEntry = [];
                    targetExit = [];


                    % book keeping & calculate focus index
                    stfElement.numOfBixelsPerRay(j) = numel([stfElement.ray(j).energy]);
                    currentMinimumFWHM = matRad_interp1(this.machine.meta.LUT_bxWidthminFWHM(1,:)',...
                        this.machine.meta.LUT_bxWidthminFWHM(2,:)',...
                        this.bixelWidth, ...
                        this.machine.meta.LUT_bxWidthminFWHM(2,end));
                    focusIx  =  ones(stfElement.numOfBixelsPerRay(j),1);
                    [~, vEnergyIx] = min(abs(bsxfun(@minus,[this.machine.data.energy]',...
                        repmat(stfElement.ray(j).energy,length([this.machine.data]),1))));

                    % get for each spot the focus index
                    for k = 1:stfElement.numOfBixelsPerRay(j)
                        focusIx(k) = find(this.machine.data(vEnergyIx(k)).initFocus.SisFWHMAtIso > currentMinimumFWHM,1,'first');
                    end

                    stfElement.ray(j).focusIx = focusIx';

                    %Get machine bounds
                    numParticlesPerMU = 1e6*ones(1,stfElement.numOfBixelsPerRay(j));
                    minMU = zeros(1,stfElement.numOfBixelsPerRay(j));
                    maxMU = Inf(1,stfElement.numOfBixelsPerRay(j));
                    for k = 1:stfElement.numOfBixelsPerRay(j)
                        if isfield(this.machine.data(vEnergyIx(k)),'MUdata')
                            MUdata = this.machine.data(vEnergyIx(k)).MUdata;
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

                    stfElement.ray(j).numParticlesPerMU = numParticlesPerMU;
                    stfElement.ray(j).minMU = minMU;
                    stfElement.ray(j).maxMU = maxMU;

                else % target not hit
                    stfElement.ray(j)               = [];
                    stfElement.numOfBixelsPerRay(j) = [];
                end

            end
        end

        function  postProc = initializePostProcessing(this,postProc)

            % get minimum energy per field
            minEnergy = min([postProc.ray.energy]);
            maxEnergy = max([postProc.ray.energy]);

            % get corresponding peak position
            minPeakPos  = this.machine.data(minEnergy == this.availableEnergies).peakPos;
            maxPeakPos  = this.machine.data(maxEnergy == this.availableEnergies).peakPos;

            % find set of energyies with adequate spacing


            postProc.longitudinalSpotSpacing = this.longitudinalSpotSpacing;

            tolerance              = this.longitudinalSpotSpacing/10;

            useEnergyBool = this.availablePeakPos >= minPeakPos & this.availablePeakPos <= maxPeakPos;

            ixCurr = find(useEnergyBool,1,'first');
            ixRun  = ixCurr + 1;
            ixEnd  = find(useEnergyBool,1,'last');

            while ixRun <= ixEnd
                if abs(this.availablePeakPos(ixRun)-this.availablePeakPos(ixCurr)) < ...
                        this.longitudinalSpotSpacing - tolerance
                    useEnergyBool(ixRun) = 0;
                else
                    ixCurr = ixRun;
                end
                ixRun = ixRun + 1;
            end

            for j = postProc.numOfRays:-1:1
                for k = postProc.numOfBixelsPerRay(j):-1:1
                    maskEnergy = postProc.ray(j).energy(k) == this.availableEnergies;
                    if ~useEnergyBool(maskEnergy)
                        postProc.ray(j).energy(k)         = [];
                        postProc.ray(j).focusIx(k)        = [];
                        postProc.ray(j).rangeShifter(k)   = [];
                        postProc.numOfBixelsPerRay(j) = postProc.numOfBixelsPerRay(j) - 1;
                    end
                end
                if isempty(postProc.ray(j).energy)
                    postProc.ray(j) = [];
                    postProc.numOfBixelsPerRay(j) = [];
                    postProc.numOfRays = postProc.numOfRays - 1;
                end
            end

        end
    end
end
