function Phi_an = get2DAnisotropyFunction(obj,r,theta)
% this function returns the value of the 2D anisotropy function either by
% returning the given value (from source specifiations), interpolation or
% extrapolation. Inter- and extrapolation are implemented according to the
% Rivard et al. (2007): Supplement to AAPM TG-43 update, page 2191f.
%
% Output parameters:
% Phi_an:   The 2D anisotropy function F(r, theta)
%
% Input paramaters:
% obj:      object of class Source.m
% r:        Vector containing the (radial) distance from the source center to
%           P(r,theta), with units of cm.
% theta:    Vector of polar angles between the longitudinal axis of the source
%           and the ray from the ac- tive source center to the calculation point

prec =1e-5;


Phi_an = zeros(size(r));


%intertype: 1 for r outside of tabulated range, 0 otherwise
indInside      = min(obj.AnisotropyRadialDistances)<=r & r<= max(obj.AnisotropyRadialDistances);
indInsideR     = find(indInside);
indOutsideR    = find(~indInside);


% inside: linear-linear interpolation for r inside the range and theta
ir1 = zeros(length(indInsideR),1);
ir2 = ir1;
it1 = ir1;
it2 = ir1;
for ir=1:length(indInsideR)
    ir1(ir) = find(obj.AnisotropyRadialDistances<=r(indInsideR(ir)), 1, 'last' );
    ir2(ir) = find(obj.AnisotropyRadialDistances>=r(indInsideR(ir)), 1, 'first' );

    if isempty(find(obj.AnisotropyPolarAngles<=theta(indInsideR(ir)),1, 'last' ))
        it1(ir) = 1;
    else
        it1(ir) = find(obj.AnisotropyPolarAngles<=theta(indInsideR(ir)),1, 'last' );
    end
    it2(ir) = find(obj.AnisotropyPolarAngles>=theta(indInsideR(ir)),1, 'first' );

end

%case 1
sameTdiffR = it1==it2 & ir1~=ir2;
if any(sameTdiffR)
    Q11 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(sameTdiffR),ir1(sameTdiffR)));
    Q12 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(sameTdiffR),ir2(sameTdiffR)));
    
    Phi_an(indInsideR(sameTdiffR)) = Q11+...
        ((Q12-Q11)./...
       (obj.AnisotropyRadialDistances(ir2(sameTdiffR))-...
        obj.AnisotropyRadialDistances(ir1(sameTdiffR)))').*(r(indInsideR(sameTdiffR))-...
        obj.AnisotropyRadialDistances(ir1(sameTdiffR))');
end
%case 2
sameRdiffT = ir1==ir2 & it1~=it2;
if any(sameRdiffT)
    Q11 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(sameRdiffT),ir1(sameRdiffT)));
    Q12 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it2(sameRdiffT),ir1(sameRdiffT)));
    
    Phi_an(indInsideR(sameRdiffT))= Q11+((Q12-Q11)./...
       (obj.AnisotropyPolarAngles(it2(sameRdiffT))-...
        obj.AnisotropyPolarAngles(it1(sameRdiffT)))).*(theta(indInsideR(sameRdiffT))-...
        obj.AnisotropyPolarAngles(it1(sameRdiffT)));
end
%case 3
diffRdiffT = ir1~=ir2 & it1~=it2;
if any(diffRdiffT)
    Q11 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(diffRdiffT),ir1(diffRdiffT)));
    Q12 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it2(diffRdiffT),ir1(diffRdiffT)));
    Q21 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(diffRdiffT),ir2(diffRdiffT)));
    Q22 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it2(diffRdiffT),ir2(diffRdiffT)));
    
    
    Phi_an(indInsideR(diffRdiffT)) = (Q11.*(obj.AnisotropyRadialDistances(ir2(diffRdiffT))'-r(indInsideR(diffRdiffT))).*(obj.AnisotropyPolarAngles(it2(diffRdiffT))-theta(indInsideR(diffRdiffT)))+...
        Q21.*(r(indInsideR(diffRdiffT))-obj.AnisotropyRadialDistances(ir1(diffRdiffT))').*(obj.AnisotropyPolarAngles(it2(diffRdiffT))-theta(indInsideR(diffRdiffT)))+...
        Q12.*(obj.AnisotropyRadialDistances(ir2(diffRdiffT))'-r(indInsideR(diffRdiffT))).*(theta(indInsideR(diffRdiffT))-obj.AnisotropyPolarAngles(it1(diffRdiffT)))+...
        Q22.*(r(indInsideR(diffRdiffT))-obj.AnisotropyRadialDistances(ir1(diffRdiffT))').*(theta(indInsideR(diffRdiffT))-obj.AnisotropyPolarAngles(it1(diffRdiffT))))./...
        ((obj.AnisotropyRadialDistances(ir2(diffRdiffT))-obj.AnisotropyRadialDistances(ir1(diffRdiffT)))'.*...
        (obj.AnisotropyPolarAngles(it2(diffRdiffT))-obj.AnisotropyPolarAngles(it1(diffRdiffT))));
end
%case 4
remaining = ~(diffRdiffT | sameRdiffT | sameTdiffR);
if any(remaining)
    Phi_an(indInsideR(remaining)) = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(remaining),ir1(remaining)));
end

% r < rmin and r > rmax nearest-neighbor or zeroth-order approach
ir1 = zeros(length(indOutsideR),1);
it1 = ir1;
it2 = ir1;
for ir=1:length(indOutsideR)
    if isempty(find(obj.AnisotropyPolarAngles<=theta(indOutsideR(ir)),1, 'last' ))
        it1(ir) = 1;
    else
    it1(ir) = find(obj.AnisotropyPolarAngles<=theta(indOutsideR(ir)),1, 'last' );
    end
    it2(ir) = find(obj.AnisotropyPolarAngles>=theta(indOutsideR(ir)),1, 'first' );
end

indSmallerR      = r(indOutsideR)<min(obj.AnisotropyRadialDistances);
indGreaterR      = ~indSmallerR;
ir1(indSmallerR) = 1;
ir1(indGreaterR) = length(obj.AnisotropyRadialDistances);

indSameT    = it1==it2;
if any(indSameT)
    Phi_an(indOutsideR(indSameT)) = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(indSameT),ir1(indSameT)));
end
 
 indDiffT = ~indSameT;
if any(indDiffT)
     Q11 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it1(indDiffT),ir1(indDiffT)));
     Q12 = obj.AnisotropyMatrix(sub2ind(size(obj.AnisotropyMatrix),it2(indDiffT),ir1(indDiffT)));
     
     Phi_an(indOutsideR(indDiffT)) = Q11+...
         ((Q12-Q11)./(obj.AnisotropyPolarAngles(it2(indDiffT))-obj.AnisotropyPolarAngles(it1(indDiffT)))).*(theta(indOutsideR(indDiffT))-obj.AnisotropyPolarAngles(it1(indDiffT)));
 end


%outsideRMax   = obj.AnisotropyMatrix(:,end)*ones(1,sum(r>max(obj.AnisotropyRadialDistances)));
%Phi_an(indGreaterR) = interp1(obj.AnisotropyPolarAngles,double(obj.AnisotropyMatrix(:,end)),theta(indGreaterR),'linearl');

% r < rmin: nearest-neighbor or zeroth-order approach
%outsideRMin   = obj.AnisotropyMatrix(:,1)*ones(1,sum(r<min(obj.AnisotropyRadialDistances)));
%Phi_an(indOutsideR) = interp1(obj.AnisotropyPolarAngles,double(obj.AnisotropyMatrix(:,1)),theta(indOutsideR),'linearl');

% concatenate to 1 matrix
%Phi_an        = [outsideRMin,inside,outsideRMax];

end

