function [cstHetero] = matRad_cstHeteroAutoassign(cst)

cstHetero = cst;
% automatically assign HeterogeneityCorrection to these segmentations
lungTissue={'Lung','GTV','PTV','CTV','ITV'};  

for i = 1:length(cst)
    if contains(cst{i,2},lungTissue)
        cstHetero{i,5}.HeterogeneityCorrection = 'Lung';
    end
end

end