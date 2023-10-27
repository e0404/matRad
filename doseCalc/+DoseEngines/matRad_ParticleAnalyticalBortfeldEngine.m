classdef matRad_ParticleAnalyticalBortfeldEngine < DoseEngines.matRad_ParticlePencilBeamEngineAbstract
    % matRad_DoseEngineParticlePB:
    %   Implements an engine for particle based dose calculation
    %   For detailed information see superclass matRad_DoseEngine
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % help edit

    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        possibleRadiationModes = {'protons'}
        name = 'Analytical Bortfeld Pencil-Beam';
    end

    properties (SetAccess = public, GetAccess = public)

        %calcBortfeldBragg = true;
        modeWidth = true;               % Boolean which defines a monoenergetic (0) and gaussian (1) energy spectrum
        %modeParabCyl = true;           % Boolean which defines a handmade, slower (0) and matLab, faster (1) parabolic cylinder function

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constants of the model
        % PHYSICAL CONSTANTS
        epsilon0         = 8.854*10^(-12);                   % Vacuum dielectric constant in (C^2/(N*m^2))
        electronCharge   = 1.602*10^(-19);                   % Electron charge in (C)
        avogadroNum      = 6.022*10^23;                      % Avogadro number

        % PARAMETERS OF THE TARGET MATERIAL
        massDensity = 0.997;                            % Mass Density of the medium in (g/cm^3)
        p           = 1.77;                             % Exponent in the Bragg-Kleemann rule
        alpha       = 2.2*10^(-3);                      % Material-dependent constant in the Bragg-Kleemann rule
        beta        = 0.012;                            % Slope parameter of the linear fluence reduction
        gammaNuc   = 0.6;                              % Fraction of locally absorbed energy in nuclear interactions

        Z           = 10;                               % N of electrons per molecule (water)
        MM          = 18.01;                            % Molar mass in g/mol (water)
        %n           = DoseEngines.matRad_DoseEngineAnalyticalBragg.ro*DoseEngines.matRad_DoseEngineAnalyticalBragg.NA/DoseEngines.matRad_DoseEngineAnalyticalBragg.MM;                         % Number density of molecules per cm^3
        %alpha1      = matRad_DoseEngine.el^2*matRad_DoseEngine.n*matRad_DoseEngine.Z/(4*pi*matRad_DoseEngine.epsilon0^2)/10^8;  % Bohr's formula for d(sigmaE)^2/dz

        % PARAMETERS OF THE BEAM
        phi0         = 1;                                % Primary proton fluence
        epsilonTail  = 0.1;                              % (Small) fraction of primary fluence \phi_0 contributing to the linear "tail" of the energy spectrum
        sigmaEnergy  = 0.01;                             % sigma of the gaussian energy spectrum

        radLenght    = 36.3;

    end

    methods

        function this = matRad_ParticleAnalyticalBortfeldEngine(pln)
            % Constructor
            %
            % call
            %   engine = DoseEngines.matRad_ParticleAnalyticalPencilBeamDoseEngine(ct,stf,pln,cst)
            %
            % input
            %   pln:                        matRad plan meta information struct

            this = this@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(pln);

        end

        function dij = calcDose(this,ct,cst,stf)
            % matRad particle dose calculation wrapper
            % can be automaticly called through matRad_calcDose or
            % matRad_calcParticleDose
            %
            % call
            %   dij = this.calcDose(ct,cst,stf,pln)
            %
            % input
            %   ct:             ct cube
            %   cst:            matRad cst struct
            %   stf:            matRad steering information struct
            %
            % output
            %   dij:            matRad dij struct
            %
            % References
            %   [1] http://iopscience.iop.org/0031-9155/41/8/005
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2015 the matRad development team. 
            % 
            % This file is part of the matRad project. It is subject to the license 
            % terms in the LICENSE file found in the top-level directory of this 
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
            % of the matRad project, including this file, may be copied, modified, 
            % propagated, or distributed except according to the terms contained in the 
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            matRad_cfg =  MatRad_Config.instance();
            
            % init dose calc
            [dij,ct,cst,stf] = this.calcDoseInit(ct,cst,stf);
                        
            %Progress Counter
            counter = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:length(stf) % loop over all beams

                % init beam
                dij = this.calcDoseInitBeam(dij,ct,cst,stf,i);                    

                for j = 1:stf(i).numOfRays % loop over all rays

                    if ~isempty(stf(i).ray(j).energy)
                        
                        currRay = this.computeRayGeometry(stf(i).ray(j),dij);

                        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray

                            counter = counter + 1;
                            this.bixelsPerBeam = this.bixelsPerBeam + 1;

                            % Display progress and update text only 200 times
                            if mod(this.bixelsPerBeam,max(1,round(stf(i).totalNumOfBixels/200))) == 0
                                    matRad_progress(this.bixelsPerBeam/max(1,round(stf(i).totalNumOfBixels/200)),...
                                                    floor(stf(i).totalNumOfBixels/max(1,round(stf(i).totalNumOfBixels/200))));
                            end

                            % update waitbar only 100 times if it is not closed
                            if mod(counter,round(dij.totalNumOfBixels/100)) == 0 && ishandle(this.hWaitbar)
                                waitbar(counter/dij.totalNumOfBixels,this.hWaitbar);
                            end

                            % remember beam and bixel number
                            if ~this.calcDoseDirect
                               dij.beamNum(counter)  = i;
                               dij.rayNum(counter)   = j;
                               dij.bixelNum(counter) = k;

                               % extract MU data if present (checks for downwards compatability)
                                minMU = 0;
                                if isfield(currRay,'minMU')
                                    minMU = currRay.minMU(k);
                                end
    
                                maxMU = Inf;
                                if isfield(currRay,'maxMU')
                                    maxMU = currRay.maxMU(k);
                                end
    
                                numParticlesPerMU = 1e6;
                                if isfield(currRay,'numParticlesPerMU')
                                    numParticlesPerMU = currRay.numParticlesPerMU(k);
                                end
    
                                dij.minMU(counter,1) = minMU;
                                dij.maxMU(counter,1) = maxMU;
                                dij.numParticlesPerMU(counter,1) = numParticlesPerMU;
                            end


                            % find energy index in base data
                            energyIx = find(this.round2(currRay.energy(k),4) == this.round2([this.machine.data.energy],4));

                            % create offset vector to account for additional offsets modelled in the base data and a potential
                            % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
                            % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
                            % interpolations in this.calcParticleDoseBixel() and avoid extrapolations.
                            offsetRadDepth = this.machine.data(energyIx).offset - currRay.rangeShifter(k).eqThickness;

                            % find depth depended lateral cut off
                            if this.dosimetricLateralCutOff == 1
                                currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth;
                            elseif this.dosimetricLateralCutOff < 1 && this.dosimetricLateralCutOff > 0
                                % perform rough 2D clipping
                                currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                    currRay.radialDist_sq <= max(this.machine.data(energyIx).LatCutOff.CutOff.^2);

                                %currIx = currRay.radDepths <= this.machine.data(energyIx).depths(end) + offsetRadDepth & ...
                                %    currRay.latDistsX.^2 + currRay.latDistsZ.^2 <= max(this.machine.data(energyIx).LatCutOff.CutOff.^2);

                                % peform fine 2D clipping
                                if length(this.machine.data(energyIx).LatCutOff.CutOff) > 1
                                    currIx(currIx) = matRad_interp1((this.machine.data(energyIx).LatCutOff.depths + offsetRadDepth)',...
                                        (this.machine.data(energyIx).LatCutOff.CutOff.^2)', currRay.radDepths(currIx)) >= currRay.radialDist_sq(currIx);
                                end
                            else
                                matRad_cfg.dispError('Cutoff must be a value between 0 and 1!')
                            end

                            % empty bixels may happen during recalculation of error
                            % scenarios -> skip to next bixel
                            if ~any(currIx)
                                continue;
                            end

                            % adjust radDepth according to range shifter
                            currRadDepths = currRay.radDepths(currIx) + currRay.rangeShifter(k).eqThickness;

                            % select correct initial focus sigma squared
                            currSigmaIni_sq = currRay.sigmaIni(k)^2;

                            % consider range shifter for protons if applicable
                            if currRay.rangeShifter(k).eqThickness > 0 && strcmp(stf(i).radiationMode,'protons')

                                % compute!
                                sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy, ...
                                    currRay.rangeShifter(k), ...
                                    currRay.SSD);

                                % add to initial sigma in quadrature
                                currSigmaIni_sq = currSigmaIni_sq +  sigmaRashi^2;

                            end
                            
                            % calculate particle dose for bixel k on ray j of beam i
                            bixelDose = this.calcParticleDoseBixel(...
                                currRadDepths, ...
                                currRay.radialDist_sq(currIx), ...
                                currSigmaIni_sq, ...
                                this.machine.data(energyIx));

                            % dij sampling is not implemented for
                            % particles If we decied to implement it,
                            % it should be implemented as interface in
                            % the superclass as well

                            %{
                                if this.enableDijSampling 
                                    [currIx,bixelDose] = this.dijSampling(currIx,bixelDose,radDepths(currIx),radialDist_sq(currIx),Type,relDoseThreshold);
                                end
                            %}

                            % Save dose for every bixel in cell array
                            this.doseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);

                            if isfield(dij,'mLETDose')
                                % calculate particle LET for bixel k on ray j of beam i
                                depths = this.machine.data(energyIx).depths + this.machine.data(energyIx).offset;
                                bixelLET = matRad_interp1(depths,this.machine.data(energyIx).LET,currRadDepths);

                                % Save LET for every bixel in cell array
                                this.letDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelLET.*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            if this.calcBioDose
                                % calculate alpha and beta values for bixel k on ray j of                  
                                [bixelAlpha, bixelBeta] = matRad_calcLQParameter(...
                                    currRadDepths,...
                                    currRay.vTissueIndex_j(currIx,:),...
                                    this.machine.data(energyIx));

                                this.alphaDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1} = sparse(this.VdoseGrid(currRay.ix(currIx)),1,bixelAlpha.*bixelDose,dij.doseGrid.numOfVoxels,1);
                                this.betaDoseTmpContainer{mod(counter-1,this.numOfBixelsContainer)+1,1}  = sparse(this.VdoseGrid(currRay.ix(currIx)),1,sqrt(bixelBeta).*bixelDose,dij.doseGrid.numOfVoxels,1);
                            end

                            %  fill the dij struct each time a
                            %  bixelContainer is calculated and at the end
                            %  of the dose calculation
                            if mod(counter,this.numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels                      
                                if this.calcDoseDirect
                                    dij = this.fillDijDirect(dij,stf,i,j,k);
                                else
                                    dij = this.fillDij(dij,stf,counter);
                                end

                            end

                        end

                    end

                end
            end

            %Close Waitbar
            if ishandle(this.hWaitbar)
                delete(this.hWaitbar);
            end
        end
    end

    methods (Access = protected)

        function dose = calcParticleDoseBixel(this, radDepths, radialDist_sq, sigmaIni_sq, baseData)
            % matRad visualization of two-dimensional dose distributions
            % on ct including segmentation
            %
            % call
            %   dose = this.calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData)
            %
            % input
            %   radDepths:      radiological depths
            %   radialDist_sq:  squared radial distance in BEV from central ray
            %   sigmaIni_sq:    initial Gaussian sigma^2 of beam at patient surface
            %   baseData:       base data required for particle dose calculation
            %
            % output
            %   dose:   particle dose at specified locations as linear vector
            %
            % References
            %   [1] http://iopscience.iop.org/0031-9155/41/8/005
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2015 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % add potential offset
            depths = baseData.depths + baseData.offset;

            % convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
            conversionFactor = 1.6021766208e-02;

            % compute integrated dose and lateral sigma
            Z = this.calcAnalyticalBragg(baseData.energy, radDepths, this.modeWidth).*conversionFactor;
            X = this.calcSigmaLatMCS(radDepths, baseData.energy);
            X = [Z, X];

            %compute lateral sigma
            sigmaSq = X(:,2).^2 + sigmaIni_sq;

            % calculate dose
            dose = baseData.LatCutOff.CompFac * exp( -radialDist_sq ./ (2*sigmaSq)) .* X(:,1) ./(2*pi*sigmaSq);

            % check if we have valid dose values
            if any(isnan(dose)) || any(dose<0)
                error('Error in particle dose calculation.');
            end
        end

        function doseVector = calcAnalyticalBragg(this, primaryEnergy, depthZ, widthMod)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   call
            %     this.calcAnalyticalBragg(PrimaryEnergy, depthz, WidthMod)
            %       ===========================================================
            %       Purpose: Compute depth-dose curve i.e. the Bragg Peak
            %                in 'Bortfeld 1998' formalism.
            %
            %       Input  : primaryEnergy -- Parameter (primaryEnergy > 0, it
            %                                 is the primary energy of the beam
            %                                 )
            %                depthZ --------- Argument (depthZ > 0,
            %                                 depth in the target material).
            %                widthMod ------- Parameter
            %                                 (widthMod = 0 for a beam that
            %                                 is monoenergetic;
            %                                 WidthMod = 1 for a beam with
            %                                 initial gaussian energy spectrum)
            %       Output : doseVector ----- Depth dose curve; same size of
            %                                 depthZ
            %       ===========================================================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       This function was inspired by the paper from
            %       Thomas Bortfeld (1997) "An analytical approximation of the
            %       Bragg curve for therapeutic proton beams".
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            numberDensity   = this.massDensity*this.avogadroNum/this.MM;                                      % Number density of molecules per cm^3
            alphaPrime      = this.electronCharge^2*numberDensity*this.Z/(4*pi*this.epsilon0^2)/10^8;   % Bohr's formula for d(sigmaE)^2/dz


            %% Conversion of depth value from mm to cm
            depthZ  = depthZ./10;

            %% Compute Range and sigma

            range = this.alpha*primaryEnergy.^this.p;                                                           % Range-Energy relation, i.e. Bragg-Kleemann rule
            sigmaMonoSquared = alphaPrime*this.p^2*this.alpha^(2/this.p)*range.^(3-2/this.p)./(3-2/this.p);     % Squared Range straggling width
            sigmaMono = sqrt(sigmaMonoSquared);                                                                 % Range straggling width

            %% Compute the width of straggling, determined by widthMod

            if widthMod == 0
                wid = sigmaMono;
            elseif widthMod ==1
                energyStraggling = this.sigmaEnergy*primaryEnergy;                                                                               % Gaussian energy straggling
                sigmaTot = sqrt( sigmaMonoSquared + (energyStraggling^2) .*(this.alpha^2) .*(this.p^2) .*(primaryEnergy.^(2*this.p-2)) );   % Total straggling contribution: range + energy
                wid = sigmaTot;
            else
                error('Wrong value for WidthMod. Choose 0 for a monoenergetic beam, or choose 1 for a gaussian energy distribution');
            end

            % COEFFICIENTS IN THE BRAGG CURVE (WITHOUT STRAGGLING)
            coeffA   = this.phi0*(1-this.epsilonTail)./(this.massDensity*this.p*this.alpha^(1/this.p)*(1+this.beta*range));
            coeffA1  = coeffA;                                                % Coefficient of D1
            coeffA2  = coeffA*this.beta*(1+this.gammaNuc*this.p);            % Coefficient of D2
            coeffA3  = coeffA*this.epsilonTail*this.p./((1-this.epsilonTail)*range);              % Coefficient of Dtail
            coeffA23 = coeffA2 + coeffA3;

            %% Definition of the Depth - Dose curve without straggling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       =======================================================
            %       Purpose: compute the depth-dose curve without energy
            %                straggling hatD(z, E) (see Bortfeld 1997)
            %       Input  : depthZ --- Argument:  depth in the target
            %                                      material
            %                range  --- Parameter: range dependent on beam
            %                                      primaryEnergy
            %       Output : hatD = hatD(z, E)
            %       =======================================================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            hatD1  = coeffA1.*(range-depthZ).^(1/this.p-1);         % MCS contribution
            hatD2  = coeffA23.*(range-depthZ).^(1/this.p);          % Nuclear and MCS contribution
            hatD   = (hatD1+hatD2)  .*(depthZ<=range)...            % Total depth-dose curve without straggling
                + 0            .*(depthZ>range);

            %% Depth-Dose curve, i.e. the Bragg peak

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       ===================================================================
            %       Purpose: compute the depth-dose curve i.e. the Bragg
            %                Peak, as defined in Bortfeld 1997 (D(z) in the
            %                paper). Using the Matlab parabolic cylinder
            %                function pu(a, z) = D_(-a-1/2)(x). In the
            %                first step, the product of gauss and parabolic
            %                cylinder function is computed.
            %       Output : depthDose = D(z) ( D(z) in Bortfeld 1997 implici-
            %                                   tly depends on E and width. )
            %       ===================================================================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            functionGaussXCyl = @(a, x) exp(-x.^2/4).*pu(-a-1/2, x);    %product of gaussian
            % and "pu"

            coeffD1     = coeffA1/wid;                                                              % coefficient
            coeffD2     = coeffA23/this.p;                                                          % coefficient
            coeffD      = wid.^(1/this.p)*gamma(1/this.p)/sqrt(2*pi);                               % coefficient
            depthDose   = coeffD.*(coeffD1.*functionGaussXCyl(-1/this.p, (depthZ-range)/wid ) ...
                + coeffD2.*functionGaussXCyl(-1/this.p-1, (depthZ-range)/wid )  );

            %% OUTPUT: compute dose vector

            % Dose is computed with hatD in the plateau region, and with
            % the parabolic cylinder function in the peak region.

            isPlateau       = depthZ <  range-10*wid;
            isPeak          = depthZ >= range-10*wid & depthZ <= range+5*wid;

            dosePlateau                     = isPlateau .* hatD;
            dosePlateau(isnan(dosePlateau)) = 0;
            dosePeak                        = isPeak    .* depthDose;
            dosePeak(isnan(dosePeak))       = 0;
            doseVector                      = dosePlateau + dosePeak;
            %end

        end

        function sigmaMCS = calcSigmaLatMCS(this, depthZ, primaryEnergy)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %call
            %   this.SigmaLatMSC_H(depthz, En)
            %       ===================================================================
            %       Purpose: Compute the lateral displacement of a particle beam due to
            %                Multiple Coulomb Scattering, as function of the depth in
            %                the target material and in Highland approximation.
            %       Input  : PrimaryEnergy -- Parameter (PrimaryEnergy > 0, it is the
            %                                 primary energy of the beam)
            %                z -------------- Argument (z > 0, it is the actual
            %                                 depth in the target material)
            %       Output : displ ---------- SigmaLatMCS_H(z, E)
            %       ===================================================================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %       This function was inspired by the paper from Gottschalk et al.1992.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %% Conversion of depth value from mm to cm
            depthZ  = depthZ./10;
            %z       = depthz;

            range = this.alpha*primaryEnergy.^this.p;        % Range-Energy relation, i.e. Bragg-Kleemann rule

            sigma1      = @(z)      14.1^2 /this.radLenght * (1+1/9*log10(z./this.radLenght)).^2;
            sigma21     = @(z)   1  ./(1-2/this.p)  .*( range.^(1-2/this.p).*(range-z).^2 - (range-z).^(3-2/this.p) );
            sigma22     = @(z)   -2*(range-z)  ./(2-2/this.p)  .*( range.^(2-2/this.p) - (range-z).^(2-2/this.p) );
            sigma23     = @(z)   1   ./(3-2/this.p)  .*( range.^(3-2/this.p) - (range-z).^(3-2/this.p) );
            sigmaTot    = @(z) this.alpha^(1/this.p)/2  *sqrt( sigma1(z) .*( sigma21(z) + sigma22(z) + sigma23(z) ) );
            sigmaBeyond   = sigmaTot(range);


            isBelowR    = depthZ<=range;
            isBeyondR   = depthZ>range;
            sigmaBelowRange     = sigmaTot(depthZ)  .* isBelowR;
            sigmaBeyondRange    = sigmaBeyond       .*isBeyondR ;

            sigmaMCS = 10.* (sigmaBelowRange + sigmaBeyondRange);      %output in mm

            sigmaMCS(depthZ==0) = 0;

        end
    
        function calcLateralParticleCutOff(this,cutOffLevel,stf)
            % matRad function to calculate a depth dependend lateral cutoff
            % for each pristine particle beam
            %
            % call
            %   this.calcLateralParticleCutOff(cutOffLevel,stf)
            %
            % input
            %   this:        current engine object includes machine base data file
            %   cutOffLevel:    cut off level - number between 0 and 1
            %   stf:          	matRad steering information struct
            %
            % output
            %   machine:    	changes in the object property machine base data file including an additional field representing the lateral
            %                    cutoff
            %
            % References
            %   -
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Copyright 2015 the matRad development team.
            %
            % This file is part of the matRad project. It is subject to the license
            % terms in the LICENSE file found in the top-level directory of this
            % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
            % of the matRad project, including this file, may be copied, modified,
            % propagated, or distributed except according to the terms contained in the
            % LICENSE file.
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if cutOffLevel <= 0.98
                warning('a lateral cut off below 0.98 may result in an inaccurate dose calculation')
            end

            conversionFactor = 1.6021766208e-02;

            % function handle for calculating depth dose for APM
            sumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);

            if (cutOffLevel < 0 || cutOffLevel > 1)
                warning('lateral cutoff is out of range - using default cut off of 0.99')
                cutOffLevel = 0.99;
            end
            % define some variables needed for the cutoff calculation
            vX = [0 logspace(-1,3,1200)]; % [mm]

            % integration steps
            r_mid          = 0.5*(vX(1:end-1) +  vX(2:end))'; % [mm]
            dr             = (vX(2:end) - vX(1:end-1))';
            radialDist_sq  = r_mid.^2;

            % number of depth points for which a lateral cutoff is determined
            numDepthVal    = 35;

            % helper function for energy selection
            round2 = @(a,b)round(a*10^b)/10^b;

            % extract SSD for each bixel
            vSSD = ones(1,length([stf.ray(:).energy]));
            cnt = 1;
            for i  = 1:length(stf.ray)
                vSSD(cnt:cnt+numel([stf.ray(i).energy])-1) = stf.ray(i).SSD;
                cnt = cnt + numel(stf.ray(i).energy);
            end

            % setup energy, focus index, sigma look up table - only consider unique rows
            [energySigmaLUT,ixUnique]  = unique([[stf.ray(:).energy]; [stf.ray(:).focusIx] ; vSSD]','rows');
            rangeShifterLUT = [stf.ray(:).rangeShifter];
            rangeShifterLUT = rangeShifterLUT(1,ixUnique);

            % find the largest inital beam width considering focus index, SSD and range shifter for each individual energy
            for i = 1:size(energySigmaLUT,1)

                % find index of maximum used energy (round to keV for numerical reasons
                energyIx = max(round2(energySigmaLUT(i,1),4)) == round2([this.machine.data.energy],4);

                currFoci = energySigmaLUT(i,2);
                sigmaIni = matRad_interp1(this.machine.data(energyIx).initFocus.dist(currFoci,:)',...
                    this.machine.data(energyIx).initFocus.sigma(currFoci,:)',...
                    energySigmaLUT(i,3));
                sigmaIni_sq = sigmaIni^2;

                % consider range shifter for protons if applicable
                if  strcmp(this.machine.meta.radiationMode,'protons') && rangeShifterLUT(i).eqThickness > 0  && ~strcmp(this.machine.meta.machine,'Generic')

                    %get max range shift
                    sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy, ...
                        rangeShifterLUT(i), ...
                        energySigmaLUT(i,3));

                    % add to initial sigma in quadrature
                    sigmaIni_sq = sigmaIni_sq +  sigmaRashi.^2;

                end

                energySigmaLUT(i,4) = sigmaIni_sq;

            end

            % find for each individual energy the broadest inital beam width
            uniqueEnergies                = unique(energySigmaLUT(:,1));
            largestSigmaSq4uniqueEnergies = NaN * ones(numel(uniqueEnergies),1);
            ix_Max                        = NaN * ones(numel(uniqueEnergies),1);
            for i = 1:numel(uniqueEnergies)
                [largestSigmaSq4uniqueEnergies(i), ix_Max(i)] = max(energySigmaLUT(uniqueEnergies(i) == energySigmaLUT(:,1),4));
            end

            % get energy indices for looping
            vEnergiesIx = find(ismember([this.machine.data(:).energy],uniqueEnergies(:,1)));
            cnt         = 0;

            % loop over all entries in the machine.data struct
            for energyIx = vEnergiesIx

                % set default depth cut off - finite value will be set during first iteration
                depthDoseCutOff = inf;

                % get the current integrated depth dose profile
                if isstruct(this.machine.data(energyIx).Z)
                    idd_org = sumGauss(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd_org = this.machine.data(energyIx).Z * conversionFactor;
                end

                [~,peakIxOrg] = max(idd_org);

                % get indices for which a lateral cutoff should be calculated
                cumIntEnergy = cumtrapz(this.machine.data(energyIx).depths,idd_org);

                peakTailRelation   = 0.5;
                numDepthValToPeak  = ceil(numDepthVal*peakTailRelation);                                                                          % number of depth values from 0 to peak position
                numDepthValTail    = ceil(numDepthVal*(1-peakTailRelation));                                                                      % number of depth values behind peak position
                energyStepsToPeak  = cumIntEnergy(peakIxOrg)/numDepthValToPeak;
                energyStepsTail    = (cumIntEnergy(end)-cumIntEnergy(peakIxOrg))/numDepthValTail;
                % make sure to include 0, peak position and end position
                vEnergySteps       = unique([0:energyStepsToPeak:cumIntEnergy(peakIxOrg) cumIntEnergy(peakIxOrg) ...
                    cumIntEnergy(peakIxOrg+1):energyStepsTail:cumIntEnergy(end) cumIntEnergy(end)]);

                [cumIntEnergy,ix] = unique(cumIntEnergy);
                depthValues       = matRad_interp1(cumIntEnergy,this.machine.data(energyIx).depths(ix),vEnergySteps);

                if isstruct(this.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else%if this.calcBortfeldBragg
                    idd = this.calcAnalyticalBragg(this.machine.data(energyIx).energy,depthValues,this.modeWidth)*conversionFactor;
                    %else
                    %    idd  = matRad_interp1(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z,depthValues) * conversionFactor;
                end

                cnt = cnt +1 ;
                % % calculate dose in spot
                baseData                   = this.machine.data(energyIx);
                baseData.LatCutOff.CompFac = 1;

                for j = 1:numel(depthValues)

                    % save depth value
                    this.machine.data(energyIx).LatCutOff.depths(j) = depthValues(j);

                    if cutOffLevel == 1
                        this.machine.data(energyIx).LatCutOff.CompFac   = 1;
                        this.machine.data(energyIx).LatCutOff.CutOff(j) = Inf;
                    else

                        % calculate dose
                        dose_r = this.calcParticleDoseBixel(depthValues(j) + baseData.offset, radialDist_sq, largestSigmaSq4uniqueEnergies(cnt), baseData);

                        cumArea = cumsum(2*pi.*r_mid.*dose_r.*dr);
                        relativeTolerance = 0.5; %in [%]
                        if abs((cumArea(end)./(idd(j)))-1)*100 > relativeTolerance
                            warning('LateralParticleCutOff: shell integration is wrong !')
                        end

                        IX = find(cumArea >= idd(j) * cutOffLevel,1, 'first');
                        this.machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;

                        if isempty(IX)
                            depthDoseCutOff = Inf;
                            warning('LateralParticleCutOff: Couldnt find lateral cut off !')
                        elseif isnumeric(IX)
                            depthDoseCutOff = r_mid(IX);
                        end

                        this.machine.data(energyIx).LatCutOff.CutOff(j) = depthDoseCutOff;

                    end
                end
            end

            %% visualization
            if this.visBoolLateralCutOff

                % determine which pencil beam should be plotted
                subIx    = ceil(numel(vEnergiesIx)/2);
                energyIx = vEnergiesIx(subIx);

                baseData       = this.machine.data(energyIx);
                focusIx        = energySigmaLUT(ix_Max(subIx),2);
                maxSSD         = energySigmaLUT(ix_Max(subIx),3);
                rangeShifter   = rangeShifterLUT(ix_Max(subIx));
                TmpCompFac     = baseData.LatCutOff.CompFac;
                baseData.LatCutOff.CompFac = 1;

                % plot 3D cutoff at one specific depth on a rather sparse grid
                sStep         = 0.5;
                vLatX         = -100 : sStep : 100; % [mm]
                dimX          = numel(vLatX);
                midPos        = round(length(vLatX)/2);
                [X,Y]         = meshgrid(vLatX,vLatX);

                radDepths     = [0:sStep:this.machine.data(energyIx).depths(end)] + this.machine.data(energyIx).offset;
                radialDist_sq = (X.^2 + Y.^2);
                radialDist_sq = radialDist_sq(:);
                mDose         = zeros(dimX,dimX,numel(radDepths));
                vDoseInt      = zeros(numel(radDepths),1);

                for kk = 1:numel(radDepths)

                    % calculate initial focus sigma
                    sigmaIni = matRad_interp1(this.machine.data(energyIx).initFocus.dist(focusIx,:)', ...
                        this.machine.data(energyIx).initFocus.sigma(focusIx,:)',maxSSD);
                    sigmaIni_sq = sigmaIni^2;

                    % consider range shifter for protons if applicable
                    if rangeShifter.eqThickness > 0 && strcmp(pln.radiationMode,'protons')

                        % compute!
                        sigmaRashi = matRad_calcSigmaRashi(this.machine.data(energyIx).energy,rangeShifter,maxSSD);

                        % add to initial sigma in quadrature
                        sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

                    end

                    mDose(:,:,kk) = reshape(this.calcParticleDoseBixel(radDepths(kk), radialDist_sq, sigmaIni_sq,baseData),[dimX dimX]);

                    [~,IX]           = min(abs((this.machine.data(energyIx).LatCutOff.depths + this.machine.data(energyIx).offset) - radDepths(kk)));
                    TmpCutOff        = this.machine.data(energyIx).LatCutOff.CutOff(IX);
                    vXCut            = vX(vX<=TmpCutOff);

                    % integration steps
                    r_mid_Cut        = (0.5*(vXCut(1:end-1) +  vXCut(2:end)))'; % [mm]
                    dr_Cut           = (vXCut(2:end) - vXCut(1:end-1))';
                    radialDist_sqCut = r_mid_Cut.^2;

                    dose_r_Cut       = this.calcParticleDoseBixel(radDepths(kk), radialDist_sqCut(:), sigmaIni_sq,baseData);

                    cumAreaCut = cumsum(2*pi.*r_mid_Cut.*dose_r_Cut.*dr_Cut);

                    if ~isempty(cumAreaCut)
                        vDoseInt(kk) = cumAreaCut(end);
                    end
                end

                % obtain maximum dose
                if isstruct(this.machine.data(energyIx).Z)
                    idd = sumGauss(depthValues,this.machine.data(energyIx).Z.mean,...
                        this.machine.data(energyIx).Z.width.^2,...
                        this.machine.data(energyIx).Z.weight) * conversionFactor;
                else
                    idd  = matRad_interp1(this.machine.data(energyIx).depths,this.machine.data(energyIx).Z,depthValues) * conversionFactor;
                end

                [~,peakixDepth] = max(idd);
                dosePeakPos = this.calcParticleDoseBixel(this.machine.data(energyIx).depths(peakixDepth), 0, sigmaIni_sq, baseData);

                vLevelsDose = dosePeakPos.*[0.01 0.05 0.1 0.9];
                doseSlice   = squeeze(mDose(midPos,:,:));
                figure,set(gcf,'Color',[1 1 1]);
                subplot(311),h=imagesc(squeeze(mDose(midPos,:,:)));hold on;
                set(h,'AlphaData', .8*double(doseSlice>0));
                contour(doseSlice,vLevelsDose,'LevelListMode','manual','LineWidth',2);hold on

                ax = gca;
                ax.XTickLabelMode = 'manual';
                ax.XTickLabel     = strsplit(num2str(ax.XTick*sStep + this.machine.data(energyIx).offset),' ')';
                ax.YTickLabelMode = 'manual';
                ax.YTickLabel     = strsplit(num2str(ax.YTick*sStep + this.machine.data(energyIx).offset),' ')';

                plot(1+(this.machine.data(energyIx).LatCutOff.depths)*sStep^-1,...
                    this.machine.data(energyIx).LatCutOff.CutOff * sStep^-1 + midPos,'rx');

                legend({'isodose 1%,5%,10% 90%','calculated cutoff'}) ,colorbar,set(gca,'FontSize',12),xlabel('z [mm]'),ylabel('x [mm]');

                entry = this.machine.data(energyIx);
                if isstruct(entry.Z)
                    idd = sumGauss(entry.depths,entry.Z.mean,entry.Z.width.^2,entry.Z.weight);
                else
                    idd = this.machine.data(energyIx).Z;
                end
                subplot(312),plot(this.machine.data(energyIx).depths,idd*conversionFactor,'k','LineWidth',2),grid on,hold on
                plot(radDepths - this.machine.data(energyIx).offset,vDoseInt,'r--','LineWidth',2),hold on,
                plot(radDepths - this.machine.data(energyIx).offset,vDoseInt * TmpCompFac,'bx','LineWidth',1),hold on,
                legend({'original IDD',['cut off IDD at ' num2str(cutOffLevel) '%'],'cut off IDD with compensation'},'Location','northwest'),
                xlabel('z [mm]'),ylabel('[MeV cm^2 /(g * primary)]'),set(gca,'FontSize',12)

                totEnergy        = trapz(this.machine.data(energyIx).depths,idd*conversionFactor) ;
                totEnergyCutOff  = trapz(radDepths,vDoseInt * TmpCompFac) ;
                relDiff          =  ((totEnergy/totEnergyCutOff)-1)*100;
                title(['rel diff of integral dose ' num2str(relDiff) '%']);
                baseData.LatCutOff.CompFac = TmpCompFac;

                subplot(313),
                if isfield(this.machine.data(energyIx),'sigma1')
                    yyaxis left;
                    plot(this.machine.data(energyIx).LatCutOff.depths,this.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on
                    plot(this.machine.data(energyIx).depths,(this.machine.data(energyIx).sigma1),':','LineWidth',2),grid on,hold on,ylabel('mm')
                    yyaxis right;
                    plot(this.machine.data(energyIx).depths,(this.machine.data(energyIx).sigma2),'-.','LineWidth',2),grid on,hold on,ylabel('mm')
                    legend({'Cutoff','sigma1','sigma2'});
                else
                    yyaxis left;plot(this.machine.data(energyIx).LatCutOff.depths,this.machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on,ylabel('mm')
                    yyaxis right;subplot(313),plot(this.machine.data(energyIx).depths,this.machine.data(energyIx).sigma,'LineWidth',2),grid on,hold on
                    legend({'Cutoff','sigma'});ylabel('mm')
                end

                set(gca,'FontSize',12),xlabel('z [mm]'),  ylabel('mm')

                % plot cutoff of different energies
                figure,set(gcf,'Color',[1 1 1]);
                cnt = 1;
                for i = vEnergiesIx
                    plot(this.machine.data(i).LatCutOff.depths,this.machine.data(i).LatCutOff.CutOff,'LineWidth',1.5),hold on
                    cellLegend{cnt} = [num2str(this.machine.data(i).energy) ' MeV'];
                    cnt = cnt + 1;
                end
                grid on, grid minor,xlabel('depth in [mm]'),ylabel('lateral cutoff in [mm]')
                title(['cutoff level = ' num2str(cutOffLevel)]),
                ylim = get(gca,'Ylim');    set(gca,'Ylim',[0 ylim(2)+3]),    legend(cellLegend)
            end
        end
    end

    methods (Static)

        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information

            msg = [];
            available = false;

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');

                %check modality
                checkModality = any(strcmp(DoseEngines.matRad_ParticleAnalyticalBortfeldEngine.possibleRadiationModes, machine.meta.radiationMode));

                preCheck = checkBasic && checkModality;

                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            checkMeta = all(isfield(machine.meta,{'SAD','BAMStoIsoDist','LUT_bxWidthminFWHM','dataType'}));

            dataType = machine.meta.dataType;
            if strcmp(dataType,'singleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','sigma','offset','initFocus'}));
            elseif strcmp(dataType,'doubleGauss')
                checkData = all(isfield(machine.data,{'energy','depths','Z','peakPos','weight','sigma1','sigma2','offset','initFocus'}));
            else
                checkData = false;
            end

            available = checkMeta && checkData;
        end
    end
end

