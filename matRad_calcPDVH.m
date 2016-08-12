function [dvhPoints,volume] = matRad_calcPDVH(Q_ref,doseVec,dij,cst,varargin)

% set dose points in dvh
if ~isempty(varargin)
    dvhPoints = varargin{1};
else
    dvhPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);
end

% calculate DVH in every scenario
if length(doseVec) > 1
    % use dij scenarios
    for Scen = 1:dij.numOfScenarios
    doseInVoi   = doseVec{Scen}(cst{1,4}{1});
    dvh(Scen,:) = sum(bsxfun(@ge,doseInVoi,dvhPoints))/numel(cst{1,4}{1})*100;    
    end
elseif length(doseVec) == 1
    % create scenarios with shifts
    for Scen = 1:cst{1,5}.VOIShift.ncase
        if isequal(cst{1,5}.VOIShift.shiftType,'rounded')
            dvh(Scen,:) = sum(bsxfun(@ge,doseVec{1}(cst{1,4}{1}-cst{1,5}.VOIShift.roundedShift.idxShift(Scen)),dvhPoints))/numel(cst{1,4}{1})*100;
            
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
            
            dvh(Scen,:) = sum(bsxfun(@ge,doseVecInterp,dvhPoints))/numel(cst{1,4}{1})*100;
        end

    end
end

% calculate PDVH
if length(doseVec) > 1 & isfield(dij,'ScenProb') & length(dij.ScenProb) == length(dvh(:,1))
    error('PDVH calculation for different scenario probabilities not implemented yet')
else
    % assume equiprobable scenarios
    for j = 1:length(dvhPoints)
        VolumePointsSorted = sort(dvh(:,j),'descend');
        ix                 = max([1 ceil(Q_ref*numel(VolumePointsSorted))]);
        volume(1,j)        = VolumePointsSorted(ix);
    end
end

end