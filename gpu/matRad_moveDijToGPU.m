function dij = matRad_moveDijToGPU(dij,precision)
%MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = 'double';
end

if ~strcmp(precision,'double')
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Only double precision supported since Matlab does not support single sparse matrices!');
end

quantities = {'physicalDose','mAlphaDose','mSqrtBetaDose','mLETDose'};

for i = 1:numel(quantities)
    if isfield(dij,quantities{i})
        dij.(quantities{i}) = cellfun(@gpuArray,dij.(quantities{i}),'UniformOutput',false);
    end
end

end
