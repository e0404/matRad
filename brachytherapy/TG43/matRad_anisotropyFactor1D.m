function PhiAn = matRad_anisotropyFactor1D(r,PhiAnTab, L)
    % anisotropy function interpolates tabulated data using
    % fifth order polynomial and approximates small and large distances
    % according to Rivard et al.2007: Supplement to the 2004 update of the
    % AAPM Task Group No. 43 Report Eq. (2).
    %
    % THIS FUNCTION IS NORMALLY CALLED INSIDE matRad_getDoseRate...
    %
    % call 
    %   PhiAn = matRad_anisotropyFactor1D(r,PhiAnTab, L)
    % 
    % input
    %   r: array of radial distances in cm!
    %   PhiAnTab: tabulated consensus data of gL according to the
    %             following cell structure:
    %             PhiAnTab{1} = AnisotropyFactorRadialDistance
    %             PhiAnTab{2} = AnisotropyFactorValue
    %
    % output
    %   PhiAn: array of the same shape as r and thet containing the
    %          interpolated and extrapolated values
    
    rmin = PhiAnTab{1}(1);
    rmax = PhiAnTab{1}(end);        
    p = polyfit(PhiAnTab{1},PhiAnTab{2},5);
    PhiAn = zeros(size(r));
    PhiAn(r>=rmin & r<=rmax) = polyval(p,r(r>=rmin & r<=rmax));
    PhiAn(r>rmax) = PhiAnTab{2}(end);
    PhiAn(r<rmin) = PhiAnTab{2}(1).*(atan(L./2./r(r<rmin))./(L.*r(r<rmin)))./(atan(L./2./rmin)./(L.*rmin));     
end