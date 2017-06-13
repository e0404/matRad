function [f, df, ddf] = pezzfunc2(x)

f = 3.*x.^2 - 10.*x + 10;

if nargout > 1
    df = 6.*x - 10;
end

if nargout > 2
    ddf = 6;
end