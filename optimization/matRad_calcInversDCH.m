function D = matRad_calcInversDCH(volume,Q,doseVec,numOfScenarios,cst)

% calculate dose that corresponds to volume
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:numOfScenarios
        dosePoints(Scen) = matRad_calcInversDVH(volume,doseVec{Scen});
    end
elseif length(doseVec) == 1
    % create scenarios with shifts
    for Scen = 1:numOfScenarios
        if isequal(cst{1,5}.VOIShift.shiftType,'rounded')
            dosePoints(Scen) = matRad_calcInversDVH(volume,doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.roundedShift.idxShift(Scen)));
            
        elseif isequal(cst{1,5}.VOIShift.shiftType,'linInterp')
            % lin interpolation in x
            c00 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y0Z0(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y0Z0(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c01 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y0Z1(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y0Z1(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c10 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y1Z0(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y1Z0(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
            c11 = doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X0Y1Z1(Scen)).*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen)) +...
                  doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.linInterpShift.idxShift.X1Y1Z1(Scen)).*cst{1,5}.VOIShift.linInterpShift.idxShift.x(Scen);
             
            % lin interpolation in y  
            c0  = c00.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen))+c10.*cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen);
            c1  = c01.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen))+c11.*cst{1,5}.VOIShift.linInterpShift.idxShift.y(Scen);
                
            % lin interpolation in z
            doseVecInterp = c0.*(1-cst{1,5}.VOIShift.linInterpShift.idxShift.z(Scen))+c1.*cst{1,5}.VOIShift.linInterpShift.idxShift.z(Scen);
            
            dosePoints(Scen) = matRad_calcInversDVH(volume,doseVecInterp);
        end
    end    
end

% sort dose
dosePoints = sort(dosePoints, 'descend');

ix = max([1 ceil(Q*numel(dosePoints))]);

D = dosePoints(ix);

end