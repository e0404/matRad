function ct = matRad_moveCtToGPU(ct, precision)
% MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = [];
end

if ~isempty(precision)
    if isfield(ct, 'cubeHU')
        ct.cubeHU = cellfun(@(x) cast(x, precision), ct.cubeHU, 'UniformOutput', false);
    end
    if isfield(ct, 'cube')
        ct.cube = cellfun(@(x) cast(x, precision), ct.cube, 'UniformOutput', false);
    end
end

if isfield(ct, 'cubeHU')
    ct.cubeHU = cellfun(@gpuArray, ct.cubeHU, 'UniformOutput', false);
end

if isfield(ct, 'cube')
    ct.cube = cellfun(@gpuArray, ct.cube, 'UniformOutput', false);
end

end
