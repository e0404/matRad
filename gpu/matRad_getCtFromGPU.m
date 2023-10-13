function ct = matRad_getCtFromGPU(ct)
%MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if isfield(ct,'cubeHU')
    ct.cubeHU = cellfun(@gather,ct.cube,'UniformOutput',false);
    ct.cubeHU = cellfun(@double,ct.cube,'UniformOutput',false);
end

if isfield(ct,'cube')
    ct.cube = cellfun(@gather,ct.cube,'UniformOutput',false);
    ct.cube = cellfun(@double,ct.cube,'UniformOutput',false);
end

end

