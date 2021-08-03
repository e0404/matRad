function gL = matRad_radialDoseFunction(r,gLTab)
    % matRad_radialDoseFuncrion interpolates tabulated data using
    % fifth order polynomial and approximates small and large distances
    % according to Rivard et al.: AAPM TG-43 update, p.669, Eq. (C1).
    % 
    % THIS FUNCTION IS NORMALLY CALLED INSIDE matRad_getDoseRate...
    %
    % call
    %   matRad_radialDoseFuncrion(r,gLTab)
    %
    % input
    %   r: array of radial distances in cm!
    %   gLTab: tabulated consensus data of gL according to the
    %          following cell structure:
    %          gLTab{1} = RadialDoseDistance
    %          gLTab{2} = RadialDoseValue
    %
    % output
    %   gL: array of the same shape as r containing the interpolated
    %       and extrapolated values


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
