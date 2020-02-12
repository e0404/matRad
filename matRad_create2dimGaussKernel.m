function kernel = matRad_create2dimGaussKernel(size, sigma, resolution)

if ~(sigma == 0)
    if (mod(size, 2) == 0)
        error('error: kernel size even');
    end

    x = -(size - 1)/2 * resolution(1):resolution(1):(size - 1)/2 * resolution(1);
    y = -(size - 1)/2 * resolution(2):resolution(2):(size - 1)/2 * resolution(2);
    [X,Y] = meshgrid(x,y);

    R = sqrt(X.^2 + Y.^2); 

    radialGauss = @(r, sigma) 1./(2*pi*sigma.^2) .* exp(-0.5 .* (r./sigma).^2);
    kernel = radialGauss(R, sigma);
    kernel = kernel/sum(kernel(:));
else
    kernel = zeros(size, size);
    kernel(floor(size/2)+1,floor(size/2)+1) = 1;
end

    
