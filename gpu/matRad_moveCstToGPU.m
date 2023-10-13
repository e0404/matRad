function cst = matRad_moveCstToGPU(cst,precision)
%MATRAD_MOVECTTOGPU Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    precision = 'double';
end

for i = 1:size(cst,1)

    if strcmp(precision,'single')
        cst{i,4} = cellfun(@single,cst{i,4},'UniformOutput',false);
    end
    cst{i,4} = cellfun(@gpuArray,cst{i,4},'UniformOutput',false);

end

end

