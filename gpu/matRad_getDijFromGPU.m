function dij = matRad_getDijFromGPU(dij)
%MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

quantities = {'physicalDose','mAlphaDose','mSqrtBetaDose','mLETDose'};

for i = 1:numel(quantities)
    if isfield(dij,quantities{i})
        dij.(quantities{i}) = cellfun(@gather,dij.(quantities{i}),'UniformOutput',false);
    end
end

end
