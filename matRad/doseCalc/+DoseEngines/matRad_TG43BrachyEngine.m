classdef matRad_TG43BrachyEngine < DoseEngines.matRad_DoseEngineBase
% Engine for dose calculation based on the revised AAPM protocol for
% brachytherapy Dose Calculations
%
%
% for more informations see superclass
% DoseEngines.matRad_DoseEngineBase
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        possibleRadiationModes = {'brachy'};
        name = 'TG43';
        shortName = 'TG43';
    end

    properties
        DistanceCutoff;
        TG43approximation;
    end

    methods

        function this = matRad_TG43BrachyEngine(pln)
            if nargin < 1
                pln = [];
            end

            % call superclass constructor
            this = this@DoseEngines.matRad_DoseEngineBase(pln);
        end

        function setDefaults(this)
            setDefaults@DoseEngines.matRad_DoseEngineBase(this);

            this.DistanceCutoff = 130;

            this.doseGrid.resolution.x = 5;
            this.doseGrid.resolution.y = 5;
            this.doseGrid.resolution.z = 5;
        end
    end

    methods (Access = protected)

        function dij = initDoseCalc(this,ct,cst,stf)

            for i = 1:numel(stf)
                stf(i).numOfRays = stf(i).numOfNeedles;
            end

            dij = initDoseCalc@DoseEngines.matRad_DoseEngineBase(this,ct,cst,stf);

            % dij meta information for brachy plan
            dij.numOfScenarios      = 1;
            dij.numOfBeams          = 1;
            dij.beamNum             = 1;
            dij.numOfNeedles        = stf.numOfNeedles;
            dij.numOfSeedsPerNeedle = stf.numOfSeedsPerNeedle;
            dij.totalNumOfSeeds     = dij.numOfNeedles*dij.numOfSeedsPerNeedle;
            dij.totalNumOfBixels    = dij.totalNumOfSeeds;
        end

        function dij = calcDose(this,ct,cst,stf)
            % Initialize dij
            dij = this.initDoseCalc(ct, cst, stf);

            matRad_cfg = MatRad_Config.instance();

            seedPoints.x = single(stf.seedPoints.x);
            seedPoints.y = single(stf.seedPoints.y);
            seedPoints.z = single(stf.seedPoints.z);

            [XGrid,YGrid,ZGrid] = meshgrid(dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            dosePoints.x = single(reshape(XGrid,1,[]));
            dosePoints.y = single(reshape(YGrid,1,[]));
            dosePoints.z = single(reshape(ZGrid,1,[]));

            matRad_cfg.dispInfo('\t computing distance transform... ');
            startTime = tic;

            DistanceMatrix = matRad_getDistanceMatrix(seedPoints,dosePoints);

            % ignore all distances > Cutoff for the following calculations to save time
            Ignore = DistanceMatrix.dist > this.DistanceCutoff;
            calcDistanceMatrix.x = DistanceMatrix.x(~Ignore);
            calcDistanceMatrix.y = DistanceMatrix.y(~Ignore);
            calcDistanceMatrix.z = DistanceMatrix.z(~Ignore);
            calcDistanceMatrix.dist = DistanceMatrix.dist(~Ignore);


            % now all fields of calcDistanceMatrix are n x 1 arrays!

            % update 
            matRad_cfg.dispInfo('done in %f s!\n',toc(startTime));
            this.progressUpdate(0.125,1);

            if ~isfield(this,'propDoseCalc') || ~isfield(this,'TG43approximation')
                this.TG43approximation = '2D';
            end

            if strcmp(this.TG43approximation,'2D')
                matRad_cfg.dispInfo('\t computing angle for TG43-2D... ');
                tmpTimer = tic;
                [ThetaMatrix,~] = this.getThetaMatrix(stf.templateNormal,calcDistanceMatrix);
                matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));
            end

            % update 
            this.progressUpdate(0.25,1);

            matRad_cfg.dispInfo('\t computing dose-rate for TG43-%s... ',this.TG43approximation);
            tmpTimer = tic;
            DoseRate = zeros(length(dosePoints.x),length(seedPoints.x));
            switch this.TG43approximation
                case '1D'
                    DoseRate(~Ignore) = ...
                        this.getDoseRate1D_poly(this.machine,calcDistanceMatrix.dist);
                case '2D'
                    DoseRate(~Ignore) = ...
                        this.getDoseRate2D_poly(this.machine,calcDistanceMatrix.dist,ThetaMatrix);
                otherwise
                    matRad_cfg.dispError('TG43 Approximation ''%s'' not known!',this.TG43approximation);
            end
            matRad_cfg.dispInfo('done in %f s!\n',toc(tmpTimer));

            dij.physicalDose = {DoseRate};

            % update waitbar, delete waitbar
            matRad_cfg.dispInfo('Brachytherapy dose calculation finished in %f s!\n',toc(startTime));
            this.progressUpdate(1,1);

            %Finalize dose calculation
            dij = this.finalizeDose(dij);
        end
    end

    methods ( Access = public )

        function dij = setCalcDose(this, ct, cst, stf)
            dij = this.calcDose(ct, cst, stf);
        end

        function DoseRate = getDoseRate1D_poly(this,machine,r_mm)
            % Calculation of radial dose Rate, interpolating using polynomes
            %   1D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update,
            %   page 639, Eq. 11:
            %
            % call
            %   DoseRate = matRad_TG43BrachyEngine.getDoseRate1D_poly(machine,r_mm)
            %
            % input
            %   machine:    TG43 information about the used seeds
            %   r:          radial distance array, given in mm!
            %
            % output
            %   DoseRate:   size(r) array of dose Rate in cGy/h
            %
            % comment on dimensions / units
            %   TG43 consensus data   cm, cGy, s
            %   matRad                mm, Gy, s
            %   output dimensions depend on the dimensions of air kerma strength
            %   Sk, normallyi in cGy*cm^2/h)

            % validate/ complete input arguments
            if ~isfield(machine.data,'AnisotropyFactorRadialDistance')
                matRad_cfg.dispError('machine is missing field "AnisotropyFactorRadialDistance"...you might be trying to apply the TG43 2D formalism for basedata measured for 1D formalism')
            end
            if ~isfield(machine.data,'AnisotropyFactorValue')
                matRad_cfg.dispError('machine is missing field "AnisotropyFactorValue"')
            end
            if ~isfield(machine.data,'lambda')
                matRad_cfg.dispError...
                    ('machine is missing field "lambda" (dose constant in water)')
            end
            if  machine.data.lambda < 0
                matRad_cfg.dispError('negative doseRate')
            end
            if min(r_mm(:)) < 0
                matRad_cfg.dispError('r contatins negative distances')
            end
            if ~isfield(machine.data,'SourceStrengthImplanted')
                this.machine.data.SourceStrengthImplanted = 1;
            end

            % arguments used during function
            % r: radius (within this function all radii are given in cm)
            r = 0.1*r_mm;

            % Sk: Air-kerma strength in  U...[1 U = 1 muGy*m^2/h) = 1 cGy*cm^2/h]
            Sk = machine.data.SourceStrengthImplanted;

            % lambda: Dose-rate constant in water (Lambda) in cGy/(h*U)
            lambda = machine.data.lambda;

            % L: length of line source in cm
            L = machine.data.ActiveSourceLength;

            % r0: reference radius in cm
            r0 = 1;

            % thet0: standard angle in degree
            theta0 = 90;

            % gLTab: Tabulated radial dose function \\ cell entry 1: radii; entry 2: values
            % radii in cm, values in units of g(r0)
            gLTab{1} = machine.data.RadialDoseDistance;
            gLTab{2} = machine.data.RadialDoseValue;

            % PhiAn: Tabulated anisotropy factor \\ cell entry 1: radii; entry 2: values
            % radii in cm, values unitless
            PhiAnTab{1} = machine.data.AnisotropyFactorRadialDistance;
            PhiAnTab{2} = machine.data.AnisotropyFactorValue;

            % 1D formalism
            % according to Rivard et al.: AAPM TG-43 update Eq. (11)
            gL = this.radialDoseFunction(r,gLTab);
            GL = this.geometryFunction(r,theta0,L);
            GL0 = this.geometryFunction(r0,theta0,L);
            PhiAn = this.anisotropyFactor1D(r,PhiAnTab, L);

            DoseRate = Sk * lambda * GL./GL0 .* gL .* PhiAn;

        end

        function DoseRate = getDoseRate2D_poly(this,machine,r_mm,theta)
            % Calculation of radial dose Rate, interpolating using polynomes
            %   2D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update,
            %   page 637, eq. 1
            %
            % call
            %   DoseRate = matRad_TG43BrachyEngine.getDoseRate2D_poly(machine,r_mm)
            %
            % input
            %   machine:    TG43 information about the used seeds
            %   r:          radial distance array, given in mm!
            %   theta:      polar angle in degree
            %
            % output
            %   DoseRate:   size(r) array of dose Rate in cGy/h
            %
            % comment on dimensions / units
            %   TG43 consensus data   cm, cGy, s
            %   matRad                mm, Gy, s
            %   output dimensions depend on the dimensions of air kerma strength
            %   Sk, normallyi in cGy*cm^2/h)

            matRad_cfg = MatRad_Config.instance();

            % validate/ complete input arguments
            if ~isfield(machine.data,'AnisotropyRadialDistances')
                matRad_cfg.dispError('machine is missing field "AnisotropyRadialDistances"...you might be trying to apply the TG43 1D formalism for basedata measured for 2D formalism')
            end
            if ~isfield(machine.data,'AnisotropyPolarAngles')
                matRad_cfg.dispError('machine is missing field "AnisotropyPolarAngles"')
            end
            if ~isfield(machine.data,'AnisotropyFunctionValue')
                matRad_cfg.dispError('machine is missing field "AnisotropyFunctionValue"')
            end
            if ~isfield(machine.data,'lambda')
                matRad_cfg.dispError('machine is missing field "lambda" (dose constant in water)')
            end
            if  machine.data.lambda < 0
                matRad_cfg.dispError('negative doseRate')
            end
            if ~isfield(machine.data,'AnisotropyRadialDistances')
                matRad_cfg.dispError('machine is missing field "AnisotropyRadialDistances"')
            end
            if ~isfield(machine.data,'AnisotropyPolarAngles')
                matRad_cfg.dispError('machine is missing field "AnisotropyPolarAngles"')
            end
            if ~isfield(machine.data,'AnisotropyFunctionValue')
                matRad_cfg.dispError('machine is missing field "AnisotropyFunctionValue"')
            end
            if min(r_mm(:)) < 0
                matRad_cfg.dispError('r contatins negative distances')
            end
            if ~isfield(machine.data,'ActiveSourceLength')
                matRad_cfg.dispError('machine is missing field "ActiveSourceLength", defining the source length')
            end
            if ~isfield(machine.data,'SourceStrengthImplanted')
                machine.data.SourceStrengthImplanted = 1;
            end

            % arguments used during function
            % r: radius (within this function all radii are given in cm)
            r = 0.1*r_mm;

            % Sk: Air-kerma strength in U...[1 U = 1 muGy*m^2/h) = 1 cGy*cm^2/h]
            Sk = machine.data.SourceStrengthImplanted;

            % lambda: Dose-rate constant in water (Lambda) in cGy/(h*U)
            lambda = machine.data.lambda;

            % L: length of line source in cm
            L = machine.data.ActiveSourceLength;

            % r0: standard radius in cm
            r0 = 1;

            % thet0: standard angle in degree
            theta0 = 90;

            % gLTab: Tabulated radial dose function \\ cell entry 1: radii; entry 2: values
            % radii in cm, values in units of g(r0)
            gLTab{1} = machine.data.RadialDoseDistance;
            gLTab{2} = machine.data.RadialDoseValue;

            % FTab: Tabulated 2D anisotropy function
            % \\ cell entry 1: radii; entry 2: angles; entry 3: values
            % radii in cm, angles in degree, values unitless
            FTab{1} = machine.data.AnisotropyRadialDistances;
            FTab{2} =machine.data.AnisotropyPolarAngles;
            FTab{3} = reshape(machine.data.AnisotropyFunctionValue, [numel(FTab{2}),numel(FTab{1})]);

            % 2D formalism
            % according to Rivard et al.: AAPM TG-43 update p. 637 eq. (1)
            % call interpolate functions and calculate formalism
            gL = this.radialDoseFunction(r,gLTab);
            GL = this.geometryFunction(r,theta,L);
            GL0 = this.geometryFunction(r0,theta0,L);
            if isfield(machine.data,'AnisotropyPolynomial')
                F = machine.data.AnisotropyPolynomial(r,theta);
            else
                F = DoseEngines.matRad_TG43BrachyEngine.anisotropyFunction2DInterp(r,theta,FTab);   % uses the Interp2 function for estimation of Anisotropy function ( Gamma(1mm,1%) pass rate 99.5%)
                %    F = matRad_TG43BrachyEngine.anisotropyFunction2D(r,theta,FTab);   % uses the 5th order polynomial for estimation of Anisotropy function
            end
            DoseRate = Sk * lambda * GL./GL0.*gL.*F;
        end

        function [ThetaMatrix,ThetaVector] = getThetaMatrix(this,templateNormal,DistanceMatrix)
            % getThetaMatrix gets (seed x dosepoint) matrix of relative polar angles
            %
            % call
            %   [ThetaMatrix,ThetaVector] = matRad_TG43BrachyEngine.getThetaMatrix(templateNormal,...
            %       DistanceMatrix)
            %   normally called within matRad_TG43BrachyEngine.getBrachyDose
            %   !!getDistanceMatrix needs to be called first!!
            %
            % input
            %   DistanceMatrix:     [dosePoint x seedPoint] struct with fields 'x','y',
            %                       'z' and total distance 'dist'
            %   templateNormal:     normal vector of template (its assumed that this is
            %                       the dir all seeds point to)
            %
            % output
            %   angle matrix:       rows: index of dosepoint
            %                       columns: index of deedpoint
            %                       entry: polar angles betreen seedpoints and
            %                       dosepoint in degrees
            %   angle vector:       column vector of angle matrix entries
            %
            % comment:
            %   The shape of the Theta matrix will be consistent with the shape of
            %   input fields.

            DistanceMatrix.dist(DistanceMatrix.dist == 0) = 1; %Avoid deviding by zero

            ThetaMatrix = acosd((templateNormal(1)*DistanceMatrix.x + templateNormal(2)*DistanceMatrix.y + templateNormal(3)*DistanceMatrix.z)./DistanceMatrix.dist);
            if nargout == 2
                ThetaVector = reshape(ThetaMatrix,[],1);
            end
        end
    end

    methods (Static)
        function PhiAn = anisotropyFactor1D(r,PhiAnTab, L)
            %   anisotropy function interpolates tabulated data
            %   using fifth order polynomial and approximates small and large distances
            %   according to Rivard et al.2007: Supplement to the 2004 update of the
            %   AAPM Task Group No. 43 Report Eq. (2).
            %   Normally called within matRad_TG43BrachyEngine.getDoseRate(...)
            %
            % call
            %   PhiAn = matRad_TG43BrachyEngine.anisotropyFactor1D(r,PhiAnTab, L)
            %
            % input
            %   r:          array of radial distances in cm!
            %   PhiAnTab:   tabulated consensus data of gL according to the following
            %               cell structure:
            %               PhiAnTab{1} = AnisotropyFactorRadialDistance
            %               PhiAnTab{2} = AnisotropyFactorValue
            %
            % output
            %   PhiAn:      array of the same shape as r and thet containing the
            %               interpolated and extrapolated values
            rmin = PhiAnTab{1}(1);
            rmax = PhiAnTab{1}(end);
            p = polyfit(PhiAnTab{1},PhiAnTab{2},5);
            PhiAn = zeros(size(r));
            PhiAn(r>=rmin & r<=rmax) = polyval(p,r(r>=rmin & r<=rmax));
            PhiAn(r>rmax) = PhiAnTab{2}(end);
            PhiAn(r<rmin) = PhiAnTab{2}(1).*(atan(L./2./r(r<rmin))./(L.*r(r<rmin)))./(atan(L./2./rmin)./(L.*rmin));
        end


        function F = anisotropyFunction2D(r,thet,FTab)
            %  anisotropy function interpolates tabulated
            %   data using fifth order polynomial and approximates small and large
            %   distances according to Rivard et al.: AAPM TG-43 update Eq. (C1).
            %   Normally called within matRad_TG43BrachyEngine.getDoseRate(...)
            %
            %   This function requires the multiPolyRegress code : https://de.mathworks.com/matlabcentral/fileexchange/34918-multivariate-polynomial-regression
            %
            % call
            %   F = matRad_TG43BrachyEngine.anisotropyFunction2D(r,thet,FTab)
            %
            % input
            %   r:      array of radial distances in cm
            %   thet:   array of azimuthal angles in °
            %   FTab:   tabulated consensus data of F according to the
            %           following cell structure:
            %           FTab{1} = AnisotropyRadialDistances
            %           FTab{2} = AnisotropyPolarAngles
            %           FTab{3} = AnisotropyFunctionValue
            %
            % output
            %   F:      array of the same shape as r and thet containing the
            %           interpolated and extrapolated values

            % prepare data for multivariate polynomial fit:
            [DataRGrid,DataThetGrid] = meshgrid(FTab{1},FTab{2});
            Data(:,1) = reshape(DataRGrid,[],1);
            Data(:,2) = reshape(DataThetGrid,[],1);
            Value     = reshape(FTab{3},[],1);
            p = MultiPolyRegress(Data,Value,5);

            % evaluate for input values
            F = p.PolynomialExpression(r,thet);

            % extrapolate for large and small values of r by taking the
            % interpolation of the maximal tabulated value at this angle
            % theta should be tabulated from 0° to 180°
            rmin = FTab{1}(1);
            rmax = FTab{1}(end);

            IndLarge = r > rmax;
            IndSmall = r < rmin;
            F(IndLarge) = p.PolynomialExpression(rmax,thet(IndLarge));
            F(IndSmall) = p.PolynomialExpression(rmin,thet(IndSmall));
        end


        function GL = geometryFunction(r,thet,L)
            %  calculates 2D geometry function
            %   according to Rivard et al.: AAPM TG-43, p 638 update Eq. (4)
            %   Normally called within matRad_TG43BrachyEngine.getDoseRate(...)
            %
            % call
            %   GL = matRad_TG43BrachyEngine.geometryFunction(r,thet,L)
            %
            % inputs
            %   r:              array of radial distances in cm!
            %   thet:           array of azimual angles in °
            %   Length:         length of radiation source in cm
            %
            % outputs
            %   GL(r,theta):    geometry function output

            % calculate solution
            if thet == 90
                beta = 2*atan(L./2./r);
                GL = beta./(L.*r);
            else
                GL = DoseEngines.matRad_TG43BrachyEngine.calcBeta(r,thet,L)./(L.*r.*sind(thet));
                GL(thet==0) = 1./(r(thet==0).^2-L^2/4);
                GL(thet==180) = 1./(r(thet==180).^2-L^2/4);
                GL(GL<0) = 0;
            end
        end

        function beta = calcBeta(r, theta,L)
            % calculate beta (see Rivard et al.: AAPM TG-43, p 637, Fig 1)
            % calculates beta from r[cm], theta [deg] and L[cm]
            % array inputs are allowed for theta

            r1 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(180 - theta)); % cos theorem
            r2 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(theta)); % cos theorem

            beta1 = asin(sind(180-theta).*L/2./r1); % sine theorem
            beta2 = asin(sind(theta).*L/2./r2); % sine theorem

            beta = beta1 + beta2;
        end

        function gL = radialDoseFunction(r,gLTab)
            %   interpolates tabulated data using
            %   fifth order polynomial and approximates small and large distances
            %   according to Rivard et al.: AAPM TG-43 update, p.669, Eq. (C1).
            %   Normally called within matRad_TG43BrachyEngine.TG43BrachyEngine.getDoseRate(...)
            %
            % call
            %   matRad_TG43BrachyEngine.TG43BrachyEngine.radialDoseFuncrion(r,gLTab)
            %
            % input
            %   r:      array of radial distances in cm!
            %   gLTab:  tabulated consensus data of gL according to the
            %           following cell structure:
            %           gLTab{1} = RadialDoseDistance
            %           gLTab{2} = RadialDoseValue
            %
            % output
            %   gL:     array of the same shape as r containing the interpolated
            %           and extrapolated values
            rmin = gLTab{1}(1);
            rmax = gLTab{1}(end);
            polyCoefficients = polyfit(gLTab{1},gLTab{2},5);
            gL = zeros(size(r));
            gL(r>=rmin & r<=rmax) = polyval(polyCoefficients,r(r>=rmin & r<=rmax));
            gL(r<rmin) = gLTab{2}(1);
            gL(r>rmax) = gLTab{2}(end) + ...
                (gLTab{2}(end)-gLTab{2}(end-1)) / (gLTab{1}(end)-...
                gLTab{1}(end-1)).*(r(r>rmax)-gLTab{1}(end));
        end

        function F = anisotropyFunction2DInterp(r,thet,FTab)
            %   anisotropy function interpolates tabulated
            %   data using interp2 ( interp technique TBD)
            %   Normally called within matRad_TG43BrachyEngine.TG43BrachyEngine.getDoseRate(...)
            %
            % call
            %   F = matRad_TG43BrachyEngine.TG43BrachyEngine.anisotropyFunction2D(r,thet,FTab)
            %
            % input
            %   r:      array of radial distances in cm
            %   thet:   array of azimuthal angles in ??
            %   FTab:   tabulated consensus data of F according to the
            %           following cell structure:
            %           FTab{1} = AnisotropyRadialDistances
            %           FTab{2} = AnisotropyPolarAngles
            %           FTab{3} = AnisotropyFunctionValue
            %
            % output
            %   F:      array of the same shape as r and thet containing the
            %           interpolated and extrapolated values
            [DataRGrid,DataThetGrid] = meshgrid(FTab{1},FTab{2});


            F = interp2(DataRGrid,DataThetGrid,FTab{3}, r, thet, 'linear',1.0);
            % extrapolate for large and small values of r by taking the
            % interpolation of the maximal tabulated value at this angle
            % theta should be tabulated from 0?? to 180??


            % rmin = FTab{1}(1);
            % rmax = FTab{1}(end);
            % 
            % IndLarge = r > rmax;
            % IndSmall = r < rmin;
            % rmaxGrid = rmax*ones(sum(IndLarge(:)),1);
            % rminGrid = rmin*ones(sum(IndSmall(:)),1);
            % F(IndLarge) = interp2(DataRGrid,DataThetGrid,FTab{3},rmaxGrid,double(thet(IndLarge)));
            % F(IndSmall) = interp2(DataRGrid,DataThetGrid,FTab{3},rminGrid,double(thet(IndSmall)));
        end

    end

end
