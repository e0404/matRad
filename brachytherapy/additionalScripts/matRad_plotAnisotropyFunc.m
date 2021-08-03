% load brachy_HDR machine from matRad\basedata first!
% execute from brachytherapy folder

% FTab: Tabulated 2D anisotropy function
FTab{1} = machine.data.AnisotropyRadialDistances;
FTab{2} = machine.data.AnisotropyPolarAngles;
FTab{3} = machine.data.AnisotropyFunctionValue;

% Plot data
[R,Thet] = meshgrid(FTab{1},FTab{2});
%[XGrid,YGrid] = pol2cart(deg2rad(Thet),R);
XGrid = Thet;
YGrid = R;
X = reshape(XGrid,[],1);
Y = reshape(YGrid,[],1);
F_Tab = reshape(FTab{3},[],1);
Data = plot3(X,Y,F_Tab,'.','MarkerSize',20,'Color','black','DisplayName','Tabulated Data');
legend();
xlabel('x[cm]')
ylabel('y[cm]')
title('5th order polynomial Fit of 2D anisotropy function')
hold on;

% Plot interpolation
rRes = 0.5; %cm
thetRes = 0.1; % deg
[r,thet] = meshgrid(0:rRes:15,0:thetRes:180);
% [x,y] = pol2cart(deg2rad(thet),r);
x = thet;
y = r;

% generate data
chris = Source(machine.data);
F = chris.get2DAnisotropyFunction(r,thet); %Dr Guthiers solution
%F = anisotropyFunction2D(r,thet,FTab); %fifth order polynomial

Fit = surf(x,y,F,'FaceAlpha',1,'EdgeColor','flat','DisplayName','Interpolation using 2D 5th order polynomial');
hold off;


    function F = anisotropyFunction2D(r,thet,FTab)
        % anisotropy function interpolates tabulated data using
        % fifth order polynomial and approximates small and large distances
        % according to Rivard et al.: AAPM TG-43 update Eq. (C1).
        
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