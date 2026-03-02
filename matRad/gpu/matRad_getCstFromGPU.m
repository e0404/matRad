function cst = matRad_getCstFromGPU(cst, precision)
% MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = [];
end

for i = 1:size(cst, 1)
    cst{i, 4} = cellfun(@gather, cst{i, 4}, 'UniformOutput', false);
    if ~isempty(precision)
        cst{i, 4} = cellfun(@(x) cast(x, precision), cst{i, 4}, 'UniformOutput', false);
    end
end

end
