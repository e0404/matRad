function DoseRate = getPointDose1D(machine,r,t)
%Calculation of radial dose Rate
%       1D dose rate formalism from Rivard et al. (2004): AAPM TG-43 update, 
%       page 639, Eq. 10:
% integrated with halflife over time t
% t input is given in days

%% input arguments
if nargin == 2
    t = Inf;    
end

%% Time Integration
halflife = machine.data.SourceIsotopeHalfLife;
TimeIntegrationFactor = 24*halflife/log(2)*(1-exp(-t*log(2)/halflife));
%% radial dose function
% from tabulated data using fifth order polynomial

polyCoefficients = polyfit(machine.data.RadialDoseDistance,machine.data.RadialDoseValue,5);
radialDoseFunction = polyval(polyCoefficients,r);

%missing: special treatment for r>r_max or r<r_min according tp AAPM report
%% 1D anisptropy function
p = polyfit(machine.data.AnisotropyFactorRadialDistance,machine.data.AnisotropyFactorValue,5);
anisotropyFunction = polyval(p,r);
% distances higher then maximal datapoint gets value of this datapoint
rMax = max(machine.data.AnisotropyFactorRadialDistance);
anisotropyFunction(r>rMax) = machine.data.AnisotropyFactorValue(end);

%% Dose Rate 1D Point source
% Sk: Air-kerma strength
Sk     = machine.data.SourceStrengthImplanted;
% lambda: Dose-rate constant in water (Lambda)
lambda = machine.data.lambda;
% r0: standard radius
r0 = 10; %mm

%dose calc within measured range
DoseRate = TimeIntegrationFactor*Sk*lambda*(r0./r).^2.*radialDoseFunction.*anisotropyFunction;
DoseRate(r>max(machine.data.RadialDoseDistance)) = 0;
end
