function cst = matRad_getCstFromGPU(cst)
%MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

for i = 1:size(cst,1)
    cst{i,4} = cellfun(@gather,cst{i,4},'UniformOutput',false);
    cst{i,4} = cellfun(@double,cst{i,4},'UniformOutput',false);
end

end

