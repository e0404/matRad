function [ddd] = matRad_getDDDfromAnalyCalc(Identifier,R0,vDepth)
% matRad function to compute the depth dose analytically according to
% T. Bortfeld's paper and code
% 
% call
%   ddd = matRad_GetDDDfromAnalyCalc(Identifier, R0, vDepth)
%
% input
%   Identifier:         string to identify particle, at the moment only
%                       'p', 'P', 'h', 'H' for protons is supported
%   R0:                 Range of the beam [in cm]
%   vDepth:             water-equivalent depths to evaluate on [in cm]
%
% output
%   ddd:                struct containing depth dose and LETd evalutaed at
%                       vDepth
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/9434986
%   [2] https://gray.mgh.harvard.edu/attachments/article/293/BraggCurve.py
%   [3] https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions
%
% Dependencies
%   Prabolic Cylinder Function Implementation [3]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch Identifier
    case {'p','P','h''H'}
           getDose = @(R0,depth) getProtonDose(R0,vDepth);
    case {'c','C'}
           getDose = @(R0,depth) getCarbonDose(R0,vDepth);
    otherwise
        error('unkown particle type')
end

ddd = getDose(R0,vDepth);

end


function [DDD] = getProtonDose(R0,depth)
% Caclulate the depth dose curve for a proton beam with given R0 and depth
% values.
% Based on T. Bortfeld code and [1]
% Edited by Mark Bangert & Niklas Wahl
% Version 0.2
p = 1.77;        % no units
%alpha = 0.0022* 3^(1-p)/2^2;  % cm/MeV^(-p)
alpha = 0.0022;
beta  = 0.012;    % 1/cm
rho   = 1.0;       % density of water g/cm^3

RT       = 0.0002;     % ??? used only for the LET computation, acts as a shift on the depth
epsilon  = 0.2;      % no tail in the initial energy spectrum
vargamma = 0.6;      % Necessary?

% Specifiy special parameters for this case

relEnergySigma = 0.014; %0.27;    % energy invariance of machine

% the original file specifies a filename at this point - not required here

%% Calculation of dose, LET_d and LET_t

% Calculate sigma
sigmaMono = 0.012 * R0^(0.935);  % ???

E0 = (R0/alpha)^(1/p);          % MeV (initial energy)
sigmaE = E0*relEnergySigma;

sigma = sqrt( sigmaMono^2+(sigmaE*alpha*p*E0^(p-1))^2); % ???

% Some more precalculations
z = depth;
zeta = (z-R0)./sigma;   % benannt wie im Paper
xi = (z-RT-R0)./sigma;  % ???

% Calculate physical DOSE
dose = sigma.^(1/p).*gamma(1/p).*(cyl_gaussj(-1./p,zeta)./sigma ...
    + cyl_gaussj(-(1/p)-1,zeta).*(beta/p+vargamma*beta+epsilon/R0) )/ ...
	(sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*R0));

	
%	calculate LETd
DDD.LETd_RT = ( 0.1 .* (sigma^(2./p).*gamma(2./p) ...
    .* ( cyl_gaussj(-2./p,xi) - cyl_gaussj(-2./p,zeta) ) ...  
    - 2 .* (RT/2.) .^ (2./p) .* exp(-(zeta+xi).^2./8)) ...
	./ (p^2.*alpha^(1./p).*(2./p-1.).*(sigma^(1./p+1.).*gamma(1./p+1) ...
    * ( cyl_gaussj(-1./p-1.,xi) -  cyl_gaussj(-1./p-1.,zeta) )  ...
    - 2.* (RT/2.).^(1./p+1.).*exp(-(zeta+xi).^2./8))));


%DDD.dose = dose./max(dose);
DDD.dose = dose;
 
end

function getCarbonDose(R0,depth)
    error('Not implemented!');
end

function y = cyl_gaussj(a,x)
%calculate product of Gaussian with parabolic cylinder function
%This is done according to [2]

branch = -12.0;

ixLargeNeg = x < branch;

%Approximation for large negative values:
x1 = x(ixLargeNeg);
y1 = sqrt(2*pi)/gamma(-a) * (-x1).^(-a-1);

y = NaN*zeros(size(x));
y(ixLargeNeg) = y1;

%
x2 = x(~ixLargeNeg);
y2a = pbdv(a,x2);
y2b = exp(-x2.^2/4);
y2 = y2a.*y2b;
y(~ixLargeNeg) = y2;


end

function val = pbdv(v,x)
    %Computes the parabolic cylinder function D_v(x) based on Weber's U(a,x)
    %Uses the package [3] https://www.mathworks.com/matlabcentral/fileexchange/22620-parabolic-cylinder-functions
    
    a = -v-0.5;
    
    ixLarge = x >= 100*abs(a);
    
    val=zeros(size(x));
    
    val(ixLarge) = arrayfun(@(x) pulx(a,x),x(ixLarge));
    val(~ixLarge) = arrayfun(@(x) pu(a,x),x(~ixLarge));
end    
    
