function ct = matRad_getCtFromGPU(ct, precision)
% MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = [];
end

if isfield(ct, 'cubeHU')
    ct.cubeHU = cellfun(@gather, ct.cube, 'UniformOutput', false);
    if ~isempty(precision)
        ct.cubeHU = cellfun(@(x) cast(x, prescision), ct.cubeHU, 'UniformOutput', false);
    end
end

if isfield(ct, 'cube')
    ct.cube = cellfun(@gather, ct.cube, 'UniformOutput', false);
    if ~isempty(precision)
        ct.cube = cellfun(@(x) cast(x, prescision), ct.cube, 'UniformOutput', false);
    end
end

end
