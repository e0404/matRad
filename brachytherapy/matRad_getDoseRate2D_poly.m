function DoseRate = matRad_getDoseRate2D_poly(machine,r_mm,thet)
% Calculation of radial dose Rate, interpolating using polynomes
%       2D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
%       page 637, eq. 1
%
% input
%   machine: TG43 information about the used seeds
%   r_mm: radial distance array, given in mm!
%   thet: polar angle in degree
%
% output 
%   DoseRate: size(r) array of dose Rate in cGy/h
%
% DIMENSIONS
% TG43 consensus data   cm, cGy/h
% matRad                mm, Gy, s

matRad_cfg = MatRad_Config.instance();

%% validate/ complete input arguments
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
if min(r_mm,[],'all') < 0
    matRad_cfg.dispError('r contatins negative distances')
end
if ~isfield(machine.data,'ActiveSourceLength')
    matRad_cfg.dispError('machine is missing field "ActiveSourceLength", defining the source length')
end
if ~isfield(machine.data,'SourceStrengthImplanted')
    machine.data.SourceStrengthImplanted = 1;
end 

%% arguments used during function
% r: radius (within this function all radii are given in cm)
r = 0.1*r_mm; 

% Sk: Air-kerma strength in U ...[1 U = 1 muGy/(m^2*h)]
Sk = machine.data.SourceStrengthImplanted;

% lambda: Dose-rate constant in water (Lambda) in cGy*/(*U)
lambda = machine.data.lambda;

% L: length of line source in cm
L = machine.data.ActiveSourceLength;

% r0: standard radius in cm
r0 = 1;

% thet0: standard angle in degree
thet0 = 90;

% gLTab: Tabulated radial dose function \\ cell entry 1: radii; entry 2: values
% radii in cm, values in units of g(r0)
gLTab{1} = machine.data.RadialDoseDistance;
gLTab{2} = machine.data.RadialDoseValue;

% FTab: Tabulated 2D anisotropy function
% \\ cell entry 1: radii; entry 2: angles; entry 3: values
% radii in cm, angles in degree, values unitless
FTab{1} = machine.data.AnisotropyRadialDistances;
FTab{2} = machine.data.AnisotropyPolarAngles;
FTab{3} = machine.data.AnisotropyFunctionValue;

% unit: factor to get from cGy/h to Gy/s
unit = 0.01/3600;

%% 2D formalism
% according to Rivard et al.: AAPM TG-43 update p. 637 eq. (1)
gL = radialDoseFuncrion(r,gLTab);
GL = geometryFunction(r,thet,L);
%GL0 = geometryFunction(r0,thet0,L); unnecessary due to built in normalization (see function)
GL0 = 1;
F = anisotropyFunction2D(r,thet,FTab);

DoseRate = unit * Sk * lambda * GL./GL0 .* gL .* F;

%% interpolation functions + geometry function
    function gL = radialDoseFuncrion(r,gLTab)
        % Radial dose function interpolates tabulated data using
        % fifth order polynomial and approximates small and large distances
        % according to Rivard et al.: AAPM TG-43 update, p.669, Eq. (C1).
        
        rmin = gLTab{1}(1);
        rmax = gLTab{1}(end);
        polyCoefficients = polyfit(gLTab{1},gLTab{2},5);
        gL = zeros(size(r));
        gL(r>=rmin & r<=rmax) = polyval(polyCoefficients,r(r>=rmin & r<=rmax));
        gL(r<rmin) = gLTab{2}(1);
        gL(r>rmax) = 0;
    end

    function GL = geometryFunction(r,thet,L) 
        % geometry function - line source approximation 
        % according to Rivard et al.: AAPM TG-43, p 638 update Eq. (4)
        %
        % INPUTS
        %   r: array of radial distances in cm
        %   Length: length of radiation source
        % OUTPUTS
        %   GL(r,theta)

        % calculate solution
        GL = calcBeta(r,thet,L)./(L.*r.*sind(thet))./... %calc according to report
                      (2*atand(L/2))*L; % normalization introduced by Claus Sarnighausen
        GL(thet==0) = 1./(r(thet==0).^2-L^2/4);
        GL(GL>10) = 10;
        GL(GL<0) = 0;
        
        
        function beta = calcBeta(r, theta,L)
            % calculate beta (see Rivard et al.: AAPM TG-43, p 637, Fig 1)
            % calculates beta from r[cm], theta [deg] and L[cm] 
            % array inputs are allowed for theta

            r1 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(180 - theta)); % cosine theorem
            r2 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(theta)); % cosine theorem

            beta1 = asind(sind(180-theta).*L/2./r1); % sine theorem
            beta2 = asind(sind(theta).*L/2./r2); % sine theorem

            beta = beta1 + beta2;
        end  
    end
        
    function F = anisotropyFunction2D(r,thet,FTab)
        % anisotropy function interpolates tabulated data using
        % fifth order polynomial and approximates small and large distances
        % according to Rivard et al.: AAPM TG-43 update Eq. (C1).
        
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

end
