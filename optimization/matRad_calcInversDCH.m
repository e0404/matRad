function D = matRad_calcInversDCH(volume,Q,doseVec,numOfScenarios,cst,dij)

% calculate dose that corresponds to volume
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
        dosePoints(Scen) = matRad_calcInversDVH(volume,doseVec{Scen});
    end
elseif length(doseVec) == 1
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
        dosePoints(Scen) = matRad_calcInversDVH(volume,doseVec{1}(cst{1,4}{1}-cst{1,5}.idxShift(Scen)));
    end    
end

% sort dose
dosePoints = sort(dosePoints, 'descend');

ix = max([1 ceil(Q*numel(dosePoints))]);

D = dosePoints(ix);

end