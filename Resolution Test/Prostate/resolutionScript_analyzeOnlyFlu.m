angularResS = [0.5 1 2 4];

%dynamic, interpolated

%percentage of the volume with at least a x% error relative to the
%reference dose
%also record objective function value
percVErr1_NN = zeros(size(angularResS));
percVErr3_NN = zeros(size(angularResS));
percVErr5_NN = zeros(size(angularResS));
percVErr10_NN = zeros(size(angularResS));
percVErr1_NN_itself = zeros(size(angularResS));
percVErr3_NN_itself = zeros(size(angularResS));
percVErr5_NN_itself = zeros(size(angularResS));
percVErr10_NN_itself = zeros(size(angularResS));
fluence_NN = zeros(size(recalc.apertureInfo.beam(1).shape(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape(1).shapeMap,2), numel(angularResS));
weight_NN = zeros(size(angularResS));
obj_NN = zeros(size(angularResS));

percVErr1_NY = zeros(size(angularResS));
percVErr3_NY = zeros(size(angularResS));
percVErr5_NY = zeros(size(angularResS));
percVErr10_NY = zeros(size(angularResS));
fluence_NY = zeros(size(recalc.apertureInfo.beam(1).shape(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape(1).shapeMap,2), numel(angularResS));
weight_NY = zeros(size(angularResS));
obj_NY = zeros(size(angularResS));

percVErr1_YN = zeros(size(angularResS));
percVErr3_YN = zeros(size(angularResS));
percVErr5_YN = zeros(size(angularResS));
percVErr10_YN = zeros(size(angularResS));
fluence_YN = zeros(size(recalc.apertureInfo.beam(1).shape(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape(1).shapeMap,2), numel(angularResS));
weight_YN = zeros(size(angularResS));
obj_YN = zeros(size(angularResS));

percVErr1_YY = zeros(size(angularResS));
percVErr3_YY = zeros(size(angularResS));
percVErr5_YY = zeros(size(angularResS));
percVErr10_YY = zeros(size(angularResS));
fluence_YY = zeros(size(recalc.apertureInfo.beam(1).shape(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape(1).shapeMap,2), numel(angularResS));
weight_YY = zeros(size(angularResS));
obj_YY = zeros(size(angularResS));

percVErr1_YY_oldDij = zeros(size(angularResS));
percVErr3_YY_oldDij = zeros(size(angularResS));
percVErr5_YY_oldDij = zeros(size(angularResS));
percVErr10_YY_oldDij = zeros(size(angularResS));
percVErr1_YY_oldDij_itself = zeros(size(angularResS));
percVErr3_YY_oldDij_itself = zeros(size(angularResS));
percVErr5_YY_oldDij_itself = zeros(size(angularResS));
percVErr10_YY_oldDij_itself = zeros(size(angularResS));
fluence_YY_oldDij = zeros(size(recalc.apertureInfo.beam(1).shape(1).shapeMap,1), size(recalc.apertureInfo.beam(1).shape(1).shapeMap,2), numel(angularResS));
weight_YY_oldDij = zeros(size(angularResS));

i = 1;
for angularRes = angularResS
    %for each angular resolution, proceed from the best approximation to
    %the worst
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('%.1f degrees, dyn + interp.mat',angularRes);
    load(fname);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_YY(:,:,i) = fluence_YY(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_YY(i) = weight_YY(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_YY(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    
    %{
    %NOT SURE IT MAKES SENSE TO DO THIS
    %next, do dynamic fluence but no interpolation
    fname = sprintf('%.1f degrees, dyn + Ninterp.mat',angularRes);
    load(fname);
    dose = recalc.resultGUI.physicalDose;
    percVErr3_YN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_YN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_YN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_YN(:,:,i) = fluence_YN(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_YN(i) = weight_YN(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_YN(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    %}
    
    %next, do interpolation but no dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + interp.mat',angularRes);
    load(fname);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_NY(:,:,i) = fluence_NY(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_NY(i) = weight_NY(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NY(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    %next, do interpolation and dynamic fluence, but using the Dij matrices
    %at the original 4degree resolution
    fname = sprintf('%.1f degrees, dyn + interp oldDij.mat',angularRes);
    load(fname);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_YY_oldDij(:,:,i) = fluence_YY_oldDij(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_YY_oldDij(i) = weight_YY_oldDij(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NY(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    %finally, do neither interpolation nor dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + Ninterp.mat',angularRes);
    load(fname);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_NN(:,:,i) = fluence_NN(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_NN(i) = weight_NN(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NN(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    
    i = i+1;
end


save('Results','*_NN', '*_NY', '*_YN', '*_YY')
