function y = matRad_interp1(xi,yi,x)

persistent isGriddedInterpolantAvailable;

if isempty(isGriddedInterpolantAvailable)
    isGriddedInterpolantAvailable = exist('griddedInterpolant','class');
end

if numel(x) == 1
    
    ix1 = find((x >= xi), 1, 'last');
    ix2 = find((x <= xi), 1, 'first');
    
    if ix1 == ix2
        y = yi(ix1);
    elseif ix2 == ix1 + 1
        y = yi(ix1) + ( yi(ix2)-yi(ix1) ) * ( x - xi(ix1) ) / ( xi(ix2) - xi(ix1) );    
    else
        error('Extrapolation not allowed');
    end
    
elseif isGriddedInterpolantAvailable
    
    if size(yi,2) > 1
        samplePoints = {xi,1:size(yi,2)};
        F = griddedInterpolant(samplePoints,yi);
    
        queryPoints = {x,1:size(yi,2)};
        y = F(queryPoints);
    else
        F = griddedInterpolant(xi,yi);
        y = F(x);
    end
else
    
    % for older matlab versions use this code
    y = interp1(xi,yi,x,'linear');

end
