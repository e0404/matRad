function [dchPoints,coverageProbabilities] = matRad_calcDCH(V_ref,doseVec,dij,cst,varargin)

% calculate inverse DVH in every scenario
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:dij.numOfScenarios
        doseInverseDVH(Scen) = matRad_calcInversDVH(V_ref,doseVec{Scen}(cst{1,4}{1}));
    end
elseif length(doseVec) == 1    
    % create scenarios with shifts
    for Scen = 1:cst{1,5}.VOIShift.ncase
        if isequal(cst{1,5}.VOIShift.shiftType,'rounded')
            doseInverseDVH(Scen) = matRad_calcInversDVH(V_ref,doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.roundedShift.idxShift(Scen)));

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
        
            doseInverseDVH(Scen) = matRad_calcInversDVH(V_ref,doseVecInterp);
        end
    end
end

% set dose points in dch, optinal: single dose point
if ~isempty(varargin)
    if isnumeric(varargin{1})
        dchPoints = varargin{1};
    elseif isequal(varargin{1},'exactDosePoints')
        dchPoints = unique(sort(doseInverseDVH));
    end
    
else
    dchPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);
end

% calculate coverage probability in every dchPoint by counting number of 
% scenarios with doseInverseDVH >= dchPoint
logicalDoseMask       = bsxfun(@ge,doseInverseDVH',dchPoints);

if length(doseVec) > 1 & isfield(dij,'ScenProb') & length(dij.ScenProb) == length(doseInverseDVH)
    % use calculated scenario probabilties
    coverageProbabilities = sum(bsxfun(@times,dij.ScenProb',logicalDoseMask))*100;
else
    % assume equiprobable scenarios
    coverageProbabilities = (1/length(doseInverseDVH))*sum(logicalDoseMask)*100;
end


end