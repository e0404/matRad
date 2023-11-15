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
        
        % target material parameters (Water)
        massDensity = 0.997;                            % Mass Density of the medium in (g/cm^3)
        p           = 1.77;                             % Exponent in the Bragg-Kleemann rule
        alpha       = 2.2*10^(-3);                      % Material-dependent constant in the Bragg-Kleemann rule
        beta        = 0.012;                            % Slope parameter of the linear fluence reduction
        gammaNuc    = 0.6;                              % Fraction of locally absorbed energy in nuclear interactions
        Z           = 10;                               % N of electrons per molecule (water)
        MM          = 18.01;                            % Molar mass in g/mol (water)

        % Beam parameters
        phi0         = 1;                               % Primary proton fluence
        epsilonTail  = 0.1;                             % (Small) fraction of primary fluence \phi_0 contributing to the linear "tail" of the energy spectrum
        sigmaEnergy  = 0.01;                            % sigma of the gaussian energy spectrum

        radLength    = 36.3;
        modeWidth   = true;                             % Boolean which defines a monoenergetic (0) and gaussian (1) energy spectrum
    end

    properties (Access = private, Constant)
        epsilon0         = 8.854*10^(-12);              % Vacuum dielectric constant in (C^2/(N*m^2))
        electronCharge   = 1.602*10^(-19);              % Electron charge in (C)
        avogadroNum      = 6.022*10^23;                 % Avogadro number
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
            
            this.calcLET = false;
        end           
    end

    methods (Access = protected)
        function [dij,ct,cst,stf] = calcDoseInit(this,ct,cst,stf)
           
            if this.calcLET == true
                matRad_cfg.dispWarning('Engine does not support LET calculation! Disabling!');
                this.calcLET = false;
            end
              
            if this.calcBioDose == true
                matRad_cfg.dispWarning('Engine does not support BioDose calculation! Disabling!');
                this.calcBioDose = false;
            end
            
            [dij,ct,cst,stf] = this.calcDoseInit@DoseEngines.matRad_ParticlePencilBeamEngineAbstract(ct,cst,stf);
        end

        function chooseLateralModel(this)
            %Now check if we need tho chose the lateral model because it
            %was set to auto
            if strcmp(this.lateralModel,'auto') 
                this.lateralModel = 'single';
            elseif ~strcmp(this.lateralModel,'single') 
                matRad_cfg.dispWarning('Engine only supports analytically computed singleGaussian lateral Model!');
                this.lateralModel = 'single';
            end              

            matRad_cfg.dispInfo('Using an analytically computed %s Gaussian pencil-beam kernel model!\n');
        end

        function [currBixel] = getBixelIndicesOnRay(this,currBixel,currRay)
            
            % create offset vector to account for additional offsets modelled in the base data and a potential
            % range shifter. In the following, we only perform dose calculation for voxels having a radiological depth
            % that is within the limits of the base data set (-> machine.data(i).dephts). By this means, we only allow
            % interpolations in this.calcParticleDoseBixel() and avoid extrapolations.
            %urrBixel.offsetRadDepth = currBixel.baseData.offset + currBixel.radDepthOffset;
            tmpOffset = currBixel.baseData.offset - currBixel.radDepthOffset;
            
            maxDepth = 1.15 * currBixel.baseData.range;

            % find depth depended lateral cut off
            if this.dosimetricLateralCutOff == 1
                currIx = currRay.radDepths <= maxDepth + tmpOffset;
            elseif this.dosimetricLateralCutOff < 1 && this.dosimetricLateralCutOff > 0
                currIx = currRay.radDepths <= maxDepth + tmpOffset;
                sigmaSq = this.calcSigmaLatMCS(currRay.radDepths(currIx) - tmpOffset, currBixel.baseData.energy).^2 + currBixel.sigmaIniSq;
                currIx(currIx) = currRay.radialDist_sq(currIx) < currBixel.baseData.LatCutOff.numSig.^2*sigmaSq;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Cutoff must be a value between 0 and 1!')
            end

            currBixel.subRayIx = currIx;
            currBixel.ix = currRay.ix(currIx);
        end

        function X = interpolateKernelsInDepth(this,bixel)
            baseData = bixel.baseData;
            
            radDepthOffset = bixel.radDepthOffset;
            
            if isfield(baseData,'offset')
                radDepthOffset = radDepthOffset - baseData.offset;
            end
           
            % calculate particle dose for bixel k on ray j of beam i
            % convert from MeV cm^2/g per primary to Gy mm^2 per 1e6 primaries
            conversionFactor = 1.6021766208e-02;
            
            energyMean = baseData.energy;
            energySpread = baseData.energy * this.sigmaEnergy;

            if isfield(baseData,'energySpectrum')
                energyMean      = baseData.energySpectrum.mean;
                energySpread    = baseData.energySpectrum.sigma/100 * baseData.energySpectrum.mean;
            end
            
            X.Z = conversionFactor * this.calcAnalyticalBragg(energyMean, bixel.radDepths + radDepthOffset, energySpread);
            X.sigma = this.calcSigmaLatMCS(bixel.radDepths, baseData.energy);
        end

        function bixel = calcParticleBixel(this,bixel)

            kernel = this.interpolateKernelsInDepth(bixel);       

            %compute lateral sigma
            sigmaSq = kernel.sigma.^2 + bixel.sigmaIniSq;

            % calculate dose
            bixel.physicalDose = bixel.baseData.LatCutOff.CompFac * exp( -bixel.radialDist_sq ./ (2*sigmaSq)) .* kernel.Z ./(2*pi*sigmaSq);

            % check if we have valid dose values
            if any(isnan(bixel.physicalDose)) || any(bixel.physicalDose<0)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Error in particle dose calculation.');
            end
        end

        function doseVector = calcAnalyticalBragg(this, primaryEnergy, depthZ, energySpread)
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
            %                energySpread ------- 
            %                                 Energy Spread
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


            % Conversion of depth value from mm to cm
            depthZ  = depthZ./10;

            % Compute Range and sigma
            range = this.alpha*primaryEnergy.^this.p;                                                           % Range-Energy relation, i.e. Bragg-Kleemann rule
            sigmaMonoSquared = alphaPrime*this.p^2*this.alpha^(2/this.p)*range.^(3-2/this.p)./(3-2/this.p);     % Squared Range straggling width
            sigmaMono = sqrt(sigmaMonoSquared);                                                                 % Range straggling width

            % Compute the width of straggling, determined by widthMod
            sigmaTot = sqrt( sigmaMonoSquared + (energySpread^2) .*(this.alpha^2) .*(this.p^2) .*(primaryEnergy.^(2*this.p-2)) );   % Total straggling contribution: range + energy
            wid = sigmaTot;


            % COEFFICIENTS IN THE BRAGG CURVE (WITHOUT STRAGGLING)
            coeffA   = this.phi0*(1-this.epsilonTail)./(this.massDensity*this.p*this.alpha^(1/this.p)*(1+this.beta*range));
            coeffA1  = coeffA;                                                % Coefficient of D1
            coeffA2  = coeffA*this.beta*(1+this.gammaNuc*this.p);            % Coefficient of D2
            coeffA3  = coeffA*this.epsilonTail*this.p./((1-this.epsilonTail)*range);              % Coefficient of Dtail
            coeffA23 = coeffA2 + coeffA3;

            % Definition of the Depth - Dose curve without straggling
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

            % Depth-Dose curve, i.e. the Bragg peak

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

            % OUTPUT: compute dose vector

            % Dose is computed with hatD in the plateau region, and with
            % the parabolic cylinder function in the peak region.

            isPlateau       = depthZ <  range-10*wid;
            isPeak          = depthZ >= range-10*wid & depthZ <= range+5*wid;

            dosePlateau                     = isPlateau .* hatD;
            dosePlateau(isnan(dosePlateau)) = 0;
            dosePeak                        = isPeak    .* depthDose;
            dosePeak(isnan(dosePeak))       = 0;
            doseVector                      = dosePlateau + dosePeak;
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


            % Conversion of depth value from mm to cm
            depthZ  = depthZ./10;
            %z       = depthz;

            range = this.alpha*primaryEnergy.^this.p;        % Range-Energy relation, i.e. Bragg-Kleemann rule

            sigma1      = @(z)      14.1^2 /this.radLength * (1+1/9*log10(z./this.radLength)).^2;
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
        
        function calcLateralParticleCutOff(this,cutOffLevel,~)
            for i = 1:numel(this.machine.data)
                this.machine.data(i).LatCutOff.CompFac = 1-cutOffLevel;
                this.machine.data(i).LatCutOff.numSig  = sqrt(2) * sqrt(gammaincinv(0.995,1)); %For a 2D symmetric gaussian we need the inverse of the incomplete Gamma function for defining the CutOff
                this.machine.data(i).LatCutOff.maxSigmaIni = max([this.machine.data(i).initFocus(:).SisFWHMAtIso]) ./ 2.3548;
                if ~isfield(this.machine.data(i),'range')
                    this.machine.data(i).range = 10 * this.alpha*this.machine.data(i).energy.^this.p;
                end
                this.machine.data(i).LatCutOff.CutOff = this.machine.data(i).LatCutOff.numSig * sqrt(this.machine.data(i).LatCutOff.maxSigmaIni^2 + this.calcSigmaLatMCS(this.machine.data(i).range,this.machine.data(i).energy)^2);
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

