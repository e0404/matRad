function DoseRate = matRad_getDoseRate1D_poly(machine,r_mm)
% Calculation of radial dose Rate, interpolating using polynomes
%       1D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
%       page 639, Eq. 11:
%
% call 
%   DoseRate = matRad_getDoseRate1D_poly(machine,r_mm)
%
% input
%   machine: TG43 information about the used seeds
%   r: radial distance array, given in mm!
%
% output 
%   DoseRate: size(r) array of dose Rate in cGy/h
%
% DIMENSIONS
% TG43 consensus data   cm, cGy, s
% matRad                mm, Gy, s

%% validate/ complete input arguments
if ~isfield(machine.data,'AnisotropyFactorRadialDistance')
    matRad_cfg.dispError('machine is missing field "AnisotropyFactorRadialDistance"...you might be trying to apply the TG43 2D formalism for basedata measured for 1D formalism') 
end
if ~isfield(machine.data,'AnisotropyFactorValue')
    matRad_cfg.dispError('machine is missing field "AnisotropyFactorValue"')
end
if ~isfield(machine.data,'lambda')
    matRad_cfg.dispError('machine is missing field "lambda" (dose constant in water)') 
end
if  machine.data.lambda < 0
    matRad_cfg.dispError('negative doseRate')
end
if min(r_mm,[],'all') < 0
    matRad_cfg.dispError('r contatins negative distances')
end
if ~isfield(machine.data,'SourceStrengthImplanted')
    machine.data.SourceStrengthImplanted = 1;
end 

%% arguments used during function
% r: radius (within this function all radii are given in cm)
r = 0.1*r_mm; 

% Sk: Air-kerma strength in muGy/(m^2*h)
Sk = machine.data.SourceStrengthImplanted;

% lambda: Dose-rate constant in water (Lambda) in cGy/(h*U)
lambda = machine.data.lambda;

% L: length of line source in cm
L = machine.data.ActiveSourceLength;

% r0: standard radius in cm
r0 = 1;

% gLTab: Tabulated radial dose function \\ cell entry 1: radii; entry 2: values
% radii in cm, values in units of g(r0)
gLTab{1} = machine.data.RadialDoseDistance;
gLTab{2} = machine.data.RadialDoseValue;

% PhiAn: Tabulated anisotropy factor \\ cell entry 1: radii; entry 2: values
% radii in cm, values unitless
PhiAnTab{1} = machine.data.AnisotropyFactorRadialDistance;
PhiAnTab{2} = machine.data.AnisotropyFactorValue;


%% 1D formalism
% according to Rivard et al.: AAPM TG-43 update Eq. (11)
gL = radialDoseFuncrion(r,gLTab);
GL = geometryFunction_thet0(r,L);
GL0 = geometryFunction_thet0(r0,L);
PhiAn = anisotropyFactor1D(r,PhiAnTab);

DoseRate = Sk * lambda * GL./GL0 .* gL .* PhiAn;

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

    function GL = geometryFunction_thet0(r,L) 
        % geometry function - line source approximation 
        % according to Rivard et al.: AAPM TG-43 update Eq. (4)
        % since in 1D we only have theta_0 = 90°, the formula simplifies to
        % GL(r,theta_0) = beta/(L*r)
        
        % INPUTS
        %   r: array of radial distances in cm
        %   Length: length of radiation source
        % OUTPUTS
        %   GL(r,theta_0) where theta_0 = 90°

        % calculate solution  
        beta = atand(L./2./r);
        GL = beta./(L.*r);  
    end
        
    function PhiAn = anisotropyFactor1D(r,PhiAnTab)
        % anisotropy function interpolates tabulated data using
        % fifth order polynomial and approximates small and large distances
        % according to Rivard et al.: AAPM TG-43 update Eq. (C1).
        rmin = PhiAnTab{1}(1);
        rmax = PhiAnTab{1}(end);        
        p = polyfit(PhiAnTab{1},PhiAnTab{2},5);
        
        PhiAn(r>=rmin & r<=rmax) = polyval(p,r(r>=rmin & r<=rmax));
        PhiAn(r>rmax) = PhiAnTab{2}(end);
        PhiAn(r<rmax) = PhiAnTab{2}(1);     
    end

end
