function dij = matRad_getDijFromGPU(dij, precision)
% MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = [];
end

quantities = {'physicalDose', 'mAlphaDose', 'mSqrtBetaDose', 'mLETDose'};

for i = 1:numel(quantities)
    if isfield(dij, quantities{i})
        dij.(quantities{i}) = cellfun(@gather, dij.(quantities{i}), 'UniformOutput', false);
        if ~isempty(precision)
            dij.(quantities{i}) = cellfun(@(x) cast(x, precision), dij.(quantities{i}), 'UniformOutput', false);
        end
    end
end

end
