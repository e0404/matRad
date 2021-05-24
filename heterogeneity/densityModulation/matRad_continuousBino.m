function [pdf,pmf,cdf] = matRad_continuousBino(input,n,p)

pdf = zeros(size(input));
pmf = pdf;
cdf = pdf;

for i = 1:size(input,1)
    X = input(i,:);
    N = n(i);
    P = p(i);
    probDens = pdf(i,:);
    propMass = pmf(i,:);
%     cumDens = cdf(i,:);
    
    % implementation according to: http://ac.inf.elte.hu/Vol_039_2013/137_39.pdf
    % for i = 1:length(X)
    % x = X(i);
    %     if x <= 0
    %         Y(i) = 0;
    %     elseif x > N+1
    %         Y(i)  = 1;
    %     elseif 0 < x <= N+1
    %         Y(i)  = betainc(x,N+1-x,p)./beta(x,N+1-x);
    %     else
    %         error('Function not defined for these input variables!')
    %     end
    % end
    
    
    % implementation according to: https://stats.stackexchange.com/questions/354214/continuous-approximation-to-binomial-distribution
    defined = (0 <= X) & (X <= N);
    continuousBino = @(k) (gamma(N+1)./gamma(k+1)./gamma(N-k+1)) .* P.^k .* (1-P).^(N-k);
    
    % calculate correction factor
    % C = ((1-p)^N-p^N)/(2*atanh(1-2*p));
    C = integral(continuousBino,0,N);
    
    continuousBinoScaled = @(k) 1/C .* continuousBino(k);
    cdfFunction = @(z) arrayfun(@(x) integral(continuousBinoScaled,0,x),z);
    
    % write distributions
    probDens(defined) = continuousBinoScaled(X(defined));
    propMass = cumsum(probDens)/sum(probDens);
%     cumDens(defined) = cdfFunction(X(defined));
%     cumDens(X>=N) = 1;
    
%     pdf(i,:) = probDens;
    pmf(i,:) = propMass;
%     cdf(i) = cumDens;

    % % Beta approximation
    % defined = (0 <= X) & (X <= 1);
    % % continuousBino = @(k) gamma(N)/gamma(p*N)/gamma((1-p)*N) * ( ((k+1)./(N+1)).^(N)/(N) - ((k)./(N+1)).^(N)/(N));
    % % Y = continuousBino(X);
    % Y(defined) = betainc(X(defined),p*N,(1-p)*N);
    % figure, plot(X,Y)
end
end
