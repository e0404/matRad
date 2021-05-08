function gL = getRadialDoseFunction(obj,r)
% This function returns the value of the radial dose function either by
% returning the given value (from source specifiations), interpolation or
% extrapolation. Inter- and extrapolation are implemented according to the
% Rivard et al. (2007): Supplement to AAPM TG-43 update. page 2195, Eq. 3.
%
% Output parameters:
% gL:       radial dose function for 1D and 2D case (???)
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to
%           P(r,theta), with units of cm.


gL = zeros(size(r));

[rMin] = min(obj.RadialDoseDistance);
[rMax] = max(obj.RadialDoseDistance);




ind     = find(rMin<=r & r<=rMax);
indMinR = zeros(length(ind),1);
indMaxR = indMinR;
for ii=1:length(ind)
    indMinR(ii) = find(obj.RadialDoseDistance<=r(ind(ii)), 1, 'last' );
    indMaxR(ii) = find(obj.RadialDoseDistance>=r(ind(ii)), 1 ,'first');
end

sameInd = indMaxR==indMinR;
diffInd = ~sameInd;

gL(ind(sameInd)) = obj.RadialDoseValue(indMinR(sameInd));
gL(ind(diffInd)) = obj.RadialDoseValue(indMinR(diffInd)).*exp(...
             (r(ind(diffInd))-obj.RadialDoseDistance(indMinR(diffInd)))./(obj.RadialDoseDistance(indMaxR(diffInd))-obj.RadialDoseDistance(indMinR(diffInd)))...
             .*(log(obj.RadialDoseValue(indMaxR(diffInd)))-log(obj.RadialDoseValue(indMinR(diffInd)))));




indSmallerR  = r<rMin;
gL(indSmallerR) = obj.RadialDoseValue(1);


indGreaterR  = r>rMax;
gL(indGreaterR) = obj.RadialDoseValue(length(obj.RadialDoseValue));





end

