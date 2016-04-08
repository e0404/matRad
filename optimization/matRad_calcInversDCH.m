function D = matRad_calcInversDCH(volume,Q,doseVec,multScen)

% calculate dose that corresponds to volume
for Scen = 1:multScen.totalNumOfScen
    dosePoints(Scen) = matRad_calcInversDVH(volume,doseVec{Scen});
end

% sort dose
dosePoints = sort(dosePoints, 'descend');

ix = max([1 ceil(Q*numel(dosePoints))]);

D = dosePoints(ix);

end