function [f, df, ddf] = pezzfunc(x)

f = 3.*x(1).^2 + 4.*x(2).^2 + x(1).*x(2) - x(1);

if nargout > 1
%     df(1) = 6.*x(1) + x(2) - 1;
%     df(2) = 8.*x(2) + x(1);
    df = [6.*x(1) + x(2) - 1, 8.*x(2) + x(1)];
    df=df';
end

if nargout > 2
%     ddf(1) = 6;
%     ddf(2) = 8;
    ddf = [6, 8];
    ddf=ddf';
end