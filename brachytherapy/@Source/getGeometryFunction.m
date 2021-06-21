function GL = getGeometryFunction(obj,r,theta)
% This function returns the value of the geometry function for the
% 2D case (Rivard et al. (2004): AAPM TG-43 update, page 638f). 
%
% Output parameters:
% GL:       geometry function 
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to 
%           P(r,theta), with units of cm.
% theta:    Vector of polar angles

%[mR,mTheta] = meshgrid(r,theta);
prec = 1e-9;

GL = zeros(size(r));

% find all places in the matrix where tetha = 0; replace by logical
% indexing
indThetaNot0 = abs(theta)>prec;
indTheta0    = ~indThetaNot0;


%% Rivard et al. (2004): AAPM TG-43 update, page 638f. Eq. 4 for theta = 0
GL(indTheta0) = (r(indTheta0).^2-obj.ActiveSourceLength^2/4).^-1;

%% Rivard et al. (2004): AAPM TG-43 update, page 638f. Eq. 4 for theta neq 0
% beta: is the angle subtended by the tips of the hypothetical line source 
% with respect to the calculation point = beta1+beta2

% beta1: angle between one tip of the line source P(r = L/2, tetha = 0°)
% and the center of the line source P(r = 0) 
% invC: 1/(distance P(r = L/2, tetha = 0°) to P(r,tetha)) from law of cosine 
% temp: from law of sines in the mentioned triangle
invC  = (sqrt((obj.ActiveSourceLength/2)^2+r(indThetaNot0).^2+2*(obj.ActiveSourceLength/2)*r(indThetaNot0).*cos(pi-theta(indThetaNot0)))).^-1;
temp  = (obj.ActiveSourceLength/2)*sin(pi-theta(indThetaNot0));
beta1 = asin(temp.*invC);

% beta2: angle between one tip of the line source P(r = L/2, tetha = 180°)
% and the center of the line source P(r = 0) 
% invA: 1/(distance P(r = L/2, tetha = 180°) to P(r,tetha)) from law of cosine 
% temp: from law of sines in the mentioned triangle
invA = (sqrt((obj.ActiveSourceLength/2)^2+r(indThetaNot0).^2+2*(obj.ActiveSourceLength/2)*r(indThetaNot0).*cos(theta(indThetaNot0)))).^-1;
temp  = (obj.ActiveSourceLength/2)*sin(theta(indThetaNot0));
beta2 = asin(temp.*invA);

% temp: part in the denomiator of the geometry function = 1/(L*r*sin(thetha))
temp  = (obj.ActiveSourceLength*r(indThetaNot0).*sin(theta(indThetaNot0))).^-1;
GL(indThetaNot0) = (beta1+beta2).*temp;

end

