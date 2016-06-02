function [dchPoints,coverageProbabilities] = matRad_calcDCH(volume,doseVec,cst,numOfScenarios,dij)

% calculate inverse DVH in every scenario
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
        doseInverseDVH(Scen) = matRad_calcInversDVH(volume/100,doseVec{Scen}(cst{1,4}{1}));
    end
elseif length(doseVec) == 1    
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
        idxShift             = cst{1,5}.voxelShift(2,Scen) + cst{1,5}.voxelShift(1,Scen)*dij.dimensions(1) + cst{1,5}.voxelShift(3,Scen)*dij.dimensions(2)*dij.dimensions(1);
        doseInverseDVH(Scen) = matRad_calcInversDVH(volume/100,doseVec{1}(cst{1,4}{1}-idxShift));
    end
end

% set dose points in dch
% dchPoints = [0,sort(doseInverseDVH),max(vertcat(doseVec{:}))*1.05];
dchPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);

% calculate coverage probability in every dchPoint by counting number of 
% scenarios with doseInverseDVH >= dchPoint
logicalDoseMask       = bsxfun(@ge,doseInverseDVH',dchPoints);
coverageProbabilities = (1/numOfScenarios)*sum(logicalDoseMask)*100;

end