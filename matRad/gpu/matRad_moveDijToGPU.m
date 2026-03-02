function dij = matRad_moveDijToGPU(dij, precision)
% MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = [];
end

quantities = {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose'};

for i = 1:numel(quantities)
    if isfield(dij, quantities{i})
        if ~isempty(precision)
            dij.(quantities{i}) = cellfun(@(x) gpuArray(cast(x, precision)), dij.(quantities{i}), 'UniformOutput', false);
        else
            dij.(quantities{i}) = cellfun(@(x) gpuArray(x), dij.(quantities{i}), 'UniformOutput', false);
        end
    end
end

end
