sigmaTot = 5;
sigmaSub = 2;
N = 5;

sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);
deltaR = (7 * sigmaHead) / N;

weight = @(mu, deltaR, i) integral(@(r) normpdf(r,mu,sigmaHead) , mu+(i-0.5)*deltaR, mu+(i+0.5)*deltaR);

range = -15:0.5:15;
sum = zeros(size(range));
allWeights = [];
for i = -(N - 1) / 2 : (N - 1) / 2
    w = weight(0, deltaR, i);
    allWeights  = [allWeights, w];
    plot(range, w * normpdf(range, i*deltaR, sigmaSub))
    sum = sum + w * normpdf(range, i*deltaR, sigmaSub);
    hold on
end
plot(range, sum)
plot(range, normpdf(range,0,sigmaTot));
hold off

radialGauss = @(r, sigma) 1 / (2 * pi * sigma^2) * exp(-0.5 * (r/sigma)^2);
[X,Y] = meshgrid(range,range);
ergGauss = arrayfun(@(x,y) radialGauss(sqrt(x^2 + y^2),sigmaTot),X,Y);
surf(ergGauss)

