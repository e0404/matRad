matRad_rc
subSigma = 2;
sigmaTot = 5;
[finalWeight, sigmaBeamlet, posX, posY, numOfSub] = matRad_calcWeights1(sigmaTot, subSigma, 21);



%%
range = -15:0.5:15;
radialGauss = @(r, sigma) 1 / (2 * pi * sigma^2) * exp(-0.5 * (r/sigma)^2);
[X,Y] = meshgrid(range,range);
ergGauss = zeros(size(X));

gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

figure
for i = 1:numOfSub
    ergGauss = ergGauss + finalWeight(i) * gauss(subSigma, X, Y, posX(i), posY(i));
        surf(X,Y,finalWeight(i) * gauss(subSigma, X, Y, posX(i), posY(i)));
        hold on
end
hold off

% realGauss = gauss(sigmaTot, X, Y, 0, 0);

% plot(ergGauss(5,:))
% hold on 
% plot(realGauss(5,:))
% surf(X,Y,ergGauss - realGauss);