sigmaTot = 5;
sigmaSub = 2;
N = 31;

sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);
deltaR = (7 * sigmaHead) / N;

weight = @(mu, deltaR, i) integral(@(r) normpdf(r,mu,sigmaHead) , mu+(i-0.5)*deltaR, mu+(i+0.5)*deltaR);

range = -15:0.1:15;
sum = zeros(size(range));
for i = -(N - 1) / 2 : (N - 1) / 2
    w = weight(0, deltaR, i);
    plot(range, w * normpdf(range, i*deltaR, sigmaSub))
    sum = sum + w * normpdf(range, i*deltaR, sigmaSub);
    hold on
end
plot(range, sum)
plot(range, normpdf(range,0,sigmaTot));
hold off