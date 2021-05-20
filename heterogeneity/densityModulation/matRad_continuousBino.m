function Y = matRad_continuousBino(X,N,p)

Y = zeros(length(X),1);
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
continuousBino = @(k) (gamma(N+1)./gamma(k((0<=k)&(k<=N))+1)./gamma(N-k((0<=k)&(k<=N))+1)) .* p.^k((0<=k)&(k<=N)) .* (1-p).^(N-k((0<=k)&(k<=N)));
% calculate correction factor
% C = ((1-p)^N-p^N)/(2*tanh(1-2*p));
C = integral(continuousBino,0,N);

Y(defined) = 1/C * continuousBino(X);

end

