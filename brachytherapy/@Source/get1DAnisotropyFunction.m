function Phi_an = get1DAnisotropyFunction(obj,r)
% this function returns the value of the 1D anisotropy function either by
% returning the given value (from source specifiations), interpolation or
% extrapolation. Inter- and extrapolation are implemented according to the 
% Rivard et al. (2007): Supplement to AAPM TG-43 update.
%
% Output parameters:
% Phi_an:   The one-dimensional anisotropy function. Dimensionless units.
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to 
%           P(r,theta), with units of cm.



Phi_an = zeros(size(r));

[rMin,indMin] = min(obj.AnisotropyFactorRadialDistance);
[rMax,indMax] = max(obj.AnisotropyFactorRadialDistance);

%if rMin==0
%    obj.AnisotropyFactorRadialDistance(indMin) = 1e-9*mm;
%    rMin = 1e-9*mm;
%    replaced = 1;
%else
%    replaced = 0;
%end
    
%%% interpolation:
% function value for distances within the range of the AnisotropyFactorRadialDistance
% log-linear approach acc. to Eq. 3 (Rivard et al. (2007): Supplement to AAPM TG-43 update, p. 2195).

%ind      = rMin<=r & r<=rMax;
%Phi_an(ind) = interp1(log(obj.AnisotropyFactorRadialDistance),obj.AnisotropyFactorValue,log(r(ind)),'linear');
% 
% ind = find(rMin<=r & r<=rMax);
% for ii=1:length(ind)
%     indMinR = find(obj.AnisotropyFactorRadialDistance<=r(ind(ii)), 1, 'last' );
%     indMaxR = find(obj.AnisotropyFactorRadialDistance>=r(ind(ii)), 1 ,'first');
%     if indMaxR==indMinR
%         Phi_an(ind(ii)) = obj.AnisotropyFactorValue(indMinR);
%     else
%         Phi_an(ind(ii)) = obj.AnisotropyFactorValue(indMinR)*exp(...
%             (r(ind(ii))-obj.AnisotropyFactorRadialDistance(indMinR))/(obj.AnisotropyFactorRadialDistance(indMaxR)-obj.AnisotropyFactorRadialDistance(indMinR))...
%             *(log(obj.AnisotropyFactorValue(indMaxR))-log(obj.AnisotropyFactorValue(indMinR))));
%     end
% end


ind     = find(rMin<=r & r<=rMax);
indMinR = zeros(length(ind),1);
indMaxR = indMinR;
for ii=1:length(ind)
    indMinR(ii) = find(obj.AnisotropyFactorRadialDistance<=r(ind(ii)), 1, 'last' );
    indMaxR(ii) = find(obj.AnisotropyFactorRadialDistance>=r(ind(ii)), 1 ,'first');
end

sameInd = indMaxR==indMinR;
diffInd = ~sameInd;

Phi_an(ind(sameInd)) = obj.AnisotropyFactorValue(indMinR(sameInd));
Phi_an(ind(diffInd)) = obj.AnisotropyFactorValue(indMinR(diffInd)).*exp(...
             (r(ind(diffInd))-obj.AnisotropyFactorRadialDistance(indMinR(diffInd)))./...
             (obj.AnisotropyFactorRadialDistance(indMaxR(diffInd))-obj.AnisotropyFactorRadialDistance(indMinR(diffInd)))...
             .*(log(obj.AnisotropyFactorValue(indMaxR(diffInd)))-log(obj.AnisotropyFactorValue(indMinR(diffInd)))));



%%% extrapolation:
% function value for distances larger than rMax, Neareast neighbor extrapolation
ind      = r>rMax;
Phi_an(ind) = obj.AnisotropyFactorValue(indMax);

% function value for distances smaller than rMin, 
% Eq. 2 (Rivard et al. (2007): Supplement to AAPM TG-43 update, p.
% 2195): Integral solved since the 1D function is independet of theta
ind      = r<rMin;
Phi_an(ind) = obj.AnisotropyFactorValue(indMin)*(atan(obj.ActiveSourceLength/2./r(ind))./atan(obj.ActiveSourceLength/2/rMin)).*(rMin./r(ind));


%if replaced
%    obj.AnisotropyFactorRadialDistance(indMin) = 0;
%end

end