function f2 = coppola(X,posx,posy,sigma_ray,x0,y0)

gaussian2 = @(x, y, mux, muy ,sig) (sqrt(2*pi)*sig)^(-1) .* exp(-(x-mux).^2/(2*(sig^2)))' * exp(-(y-muy).^2/(2*(sig^2)));
f1 = gaussian2(x0,y0,0,0,sigma_ray);
f2 = @(X,posx,posy,sigma_ray,x0,y0) -f1;
for i=1:n
    f2 = @(X,posx,posy,sigma_ray,x0,y0) f2(X,posx,posy,sigma_ray,x0,y0) + X(1) .* ...
    gaussian2(posx(i),posy(i),0,0,X(2)).*gaussian2(x0,y0,posx(i),posy(i),sigma_ray);
end