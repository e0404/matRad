 function F = matRad_anisotropyFunction2D(r,thet,FTab)
    % anisotropy function interpolates tabulated data using
    % fifth order polynomial and approximates small and large distances
    % according to Rivard et al.: AAPM TG-43 update Eq. (C1).
    %
    % THIS FUNCTION IS NORMALLY CALLED INSIDE matRad_getDoseRate...
    %
    % call
    %   F = matRad_anisotropyFunction2D(r,thet,FTab)
    %
    % input 
    %   r: array of radial distances in cm
    %   thet: array of azimuthal angles in °
    %   FTab: tabulated consensus data of F according to the
    %          following cell structure:
    %          FTab{1} = AnisotropyRadialDistances
    %          FTab{2} = AnisotropyPolarAngles
    %          FTab{3} = AnisotropyFunctionValue
    %
    % output
    %   F: array of the same shape as r and thet containing the
    %      interpolated and extrapolated values

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