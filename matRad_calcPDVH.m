function [dvhPoints,volume] = matRad_calcPDVH(coverage,doseVec,cst,numOfScenarios,dij)

% set dose points in dvh
dvhPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);

% calculate DVH in every scenario
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
    doseInVoi   = doseVec{Scen}(cst{1,4}{1});
    dvh(Scen,:) = sum(bsxfun(@ge,doseInVoi,dvhPoints))/numel(cst{1,4}{1})*100;    
    end
elseif length(doseVec) == 1
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
    idxShift  = cst{1,5}.voxelShift(2,Scen) + cst{1,5}.voxelShift(1,Scen)*dij.dimensions(1) + cst{1,5}.voxelShift(3,Scen)*dij.dimensions(2)*dij.dimensions(1);
    doseInVoi = doseVec{1}(cst{1,4}{1}-idxShift);
    dvh(Scen,:) = sum(bsxfun(@ge,doseInVoi,dvhPoints))/numel(cst{1,4}{1})*100;
    end
end

% calculate PDVH
for j = 1:length(dvhPoints)
    VolumePointsSorted = sort(dvh(:,j),'descend');
    ix                 = max([1 ceil(coverage/100*numel(VolumePointsSorted))]);
    volume(1,j)        = VolumePointsSorted(ix);
end



end