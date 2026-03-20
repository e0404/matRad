classdef matRad_ParticleSequencer < matRad_SequencerBase
    % UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        name                    = 'Particle IMPT Scanning Sequencing'
        shortName               = 'IMPT'
        possibleRadiationModes  = {'protons', 'helium', 'carbon'}
        weightPencilBeam        = 1e6
    end

    properties
        esTime = 3 * 10^6                 % [\mu s]  time required for synchrotron to recharge it' spill
        spillRechargeTime = 2 * 10^6      % [\mu s]  number of particles generated in each spill
        spillSize = 4 * 10^10
        scanSpeed = 10                    % [m/s] speed of synchrotron's lateral scanning in an IES
        spillIntensity = 4 * 10^8       % number of particles per second
    end

    methods

        function sequence = sequence(this, w, stf)

            sequence = this.calcSpotOrder(stf);
            sequence = this.calcSpotTime(sequence, w, stf);
        end

        function sequence = calcSpotOrder(~, stf)

            sequence = struct;

            wOffset = 0;
            % first loop loops over all bixels to store their position and ray number in each IES
            for i = 1:length(stf)

                usedEnergies = unique([stf(i).ray(:).energy]);
                usedEnergiesSorted = sort(usedEnergies, 'descend');

                sequence(i).orderToSTF      = zeros(stf(i).totalNumOfBixels, 1);
                sequence(i).orderToSS       = zeros(stf(i).totalNumOfBixels, 1);
                sequence(i).time            = zeros(stf(i).totalNumOfBixels, 1);
                sequence(i).e               = zeros(stf(i).totalNumOfBixels, 1);

                for e = 1:length(usedEnergies) % looping over IES's

                    s = 1;

                    for j = 1:stf(i).numOfRays % looping over all rays

                        % find the rays which are active in current IES
                        if any(stf(i).ray(j).energy == usedEnergiesSorted(e))

                            x = stf(i).ray(j).rayPos_bev(1);
                            y = stf(i).ray(j).rayPos_bev(3);

                            sequence(i).IES(e).x(s)       = x; % store x position
                            sequence(i).IES(e).y(s)       = y; % store y position
                            sequence(i).IES(e).wIndex(s) = wOffset + ...
                                sum(stf(i).numOfBixelsPerRay(1:(j - 1))) + ...
                                find(stf(i).ray(j).energy == usedEnergiesSorted(e)); % store index

                            s = s + 1;

                        end
                    end
                end

                wOffset = wOffset + sum(stf(i).numOfBixelsPerRay);

            end

        end

        function sequence = calcSpotTime(this, sequence, w, stf)
            steerTime = [stf.bixelWidth] * (10^3) / this.scanSpeed; % [\mu s]
            % after storing all the required information,
            % same loop over all bixels will put each bixel in it's order

            spillUsage = 0;
            offset = 0;

            for i = 1:length(stf)

                usedEnergies = unique([stf(i).ray(:).energy]);

                t = 0;
                orderCount = 1;

                for e = 1:length(usedEnergies)

                    % sort the y positions from high to low (backforth is up do down)
                    y_sorted = sort(unique(sequence(i).IES(e).y), 'descend');
                    x_sorted = sort(sequence(i).IES(e).x, 'ascend');

                    for k = 1:length(y_sorted)

                        y = y_sorted(k);
                        % find indexes corresponding to current y position
                        % in other words, number of bixels in the current row
                        ind_y = find(sequence(i).IES(e).y == y);

                        % since backforth fasion is zig zag like, flip the order every
                        % second row
                        if ~rem(k, 2)
                            ind_y = fliplr(ind_y);
                        end

                        % loop over all the bixels in the row
                        for is = 1:length(ind_y)

                            s = ind_y(is);

                            x = x_sorted(s);

                            wIndex = sequence(i).IES(e).wIndex(s);

                            % in case there were holes inside the plan "multi"
                            % multiplies the steertime to take it into account:
                            if k == 1 && is == 1
                                x_prev = x;
                                y_prev = y;
                            end
                            % x direction
                            multi = abs(x_prev - x) / stf(i).bixelWidth;
                            % y direction
                            multi = multi + abs(y_prev - y) / stf(i).bixelWidth;
                            %
                            x_prev = x;
                            y_prev = y;

                            % calculating the time:

                            % required spot fluence
                            numOfParticles = w(wIndex) * this.weightPencilBeam;
                            % time spent to spill the required spot fluence
                            spillTime = numOfParticles * 10^6 / this.spillIntensity;

                            % spotTime:time spent to steer scan along IES per bixel
                            t = t + multi * steerTime(i) + spillTime;

                            % taking account of the time to recharge the spill in case
                            % the required fluence was more than spill size
                            if spillUsage + numOfParticles > this.spillSize
                                t = t + this.spillRechargeTime;
                                spillUsage = 0;
                            end

                            % used amount of fluence from current spill
                            spillUsage = spillUsage + numOfParticles;

                            % storing the time and the order of bixels

                            % make the both counter and index 'per beam' - help index
                            wInd = wIndex - offset;

                            % timeline according to the spot scanning order
                            sequence(i).time(orderCount) = t;
                            % IES of bixels according to the spot scanning order
                            sequence(i).e(orderCount) = e;
                            % according to spot scanning order, sorts w index of all
                            % bixels, use this order to transfer STF order to Spot
                            % Scanning order
                            sequence(i).orderToSS(orderCount) = wInd;

                            % according to STF order, gives us order of irradiation of
                            % each bixel, use this order to transfer Spot Scanning
                            % order to STF order
                            % orderToSTF(orderToSS) = orderToSS(orderToSTF) = 1:#bixels
                            sequence(i).orderToSTF(wInd) = orderCount;

                            orderCount  = orderCount + 1;

                        end
                    end

                    t = t + this.esTime;

                end

                % storing the fluence per beam
                sequence(i).w = w(offset + 1:offset + stf(i).totalNumOfBixels);

                offset = offset + stf(i).totalNumOfBixels;
            end
        end

    end

    methods  (Static)

        function [available, msg] = isAvailable(pln, machine)
            % see superclass for information

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end
            % checkBasic
            available = isfield(machine, 'meta') && isfield(machine, 'data');

            available = available && any(isfield(machine.meta, {'machine', 'radiationMode'}));

            if ~available
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
            else
                msg = [];
            end

            % check modality
            checkModality = any(strcmp(matRad_ParticleSequencer.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_ParticleSequencer.possibleRadiationModes, pln.radiationMode));

            % Sanity check compatibility
            if checkModality
                checkModality = strcmp(machine.meta.radiationMode, pln.radiationMode);
            end

            available  = available && checkModality;

        end

        function sequence = makePhaseMatrix(sequence, numOfPhases, motionPeriod)

            phaseTime = motionPeriod * 10^6 / numOfPhases;      % time of each phase [/mu s]

            for i = 1:length(sequence)

                realTime = phaseTime;
                sequence(i).phaseMatrix = zeros(length(sequence(i).time), numOfPhases);

                iPhase = 1;
                iTime = 1;

                while iTime <= length(sequence(i).time)
                    if sequence(i).time(iTime) < realTime
                        while iTime <= length(sequence(i).time) && sequence(i).time(iTime) < realTime
                            sequence(i).phaseMatrix(iTime, iPhase) = 1;
                            iTime = iTime + 1;
                        end
                    else

                        iPhase = iPhase + 1;
                        % back to 1 after going over all phases
                        if iPhase > numOfPhases
                            iPhase = 1;
                        end
                        realTime = realTime + phaseTime;
                    end
                end

                % permuatation of phaseMatrix from SS order to STF order
                sequence(i).phaseMatrix = sequence(i).phaseMatrix(sequence(i).orderToSTF, :);
                sequence(i).phaseNum = find(sequence(i).phaseMatrix');
                % inserting the fluence in phaseMatrix
                sequence(i).phaseMatrix = sequence(i).phaseMatrix .* sequence(i).w;
            end
        end

    end

end
