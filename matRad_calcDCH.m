function [dchPoints,coverageProbabilities] = matRad_calcDCH(volume,doseVec,cst,numOfScenarios,dij)

% set dose points in DCH 
n         = 10000;
dchPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,n);

% calculate inverse DVH in every scenario and the deviation from dchPoints
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
        doseInverseDVH = matRad_calcInversDVH(volume/100,doseVec{Scen}(cst{1,4}{1}));
        dev(Scen,:)    = doseInverseDVH - dchPoints;
    end
elseif length(doseVec) == 1    
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
        idxShift       = cst{1,5}.voxelShift(2,Scen) + cst{1,5}.voxelShift(1,Scen)*dij.dimensions(1) + cst{1,5}.voxelShift(3,Scen)*dij.dimensions(2)*dij.dimensions(1);
        doseInverseDVH = matRad_calcInversDVH(volume/100,doseVec{1}(cst{1,4}{1}-idxShift));
        dev(Scen,:)    = doseInverseDVH - dchPoints;
    end
end

% calculate coverage probability in every dchPoint by counting number of 
% scenarios with doseInverseDVH >= dchPoint <-> dev >= 0
coverageProbabilities = (1/numOfScenarios)*sum(dev >= 0)*100;

end