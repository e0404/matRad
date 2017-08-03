%this is the reference plan, the most accurate way of calculating dose
fname = sprintf('0.5 degrees, dyn + interp.mat');
load(fname)
refDose = recalc.resultGUI.physicalDose;

%refDose for not interpolated, not dynamic
fname = sprintf('0.5 degrees, Ndyn + Ninterp.mat');
load(fname)
refDose_NN = recalc.resultGUI.physicalDose;

%refDose for interpolated, dynamic, old Dij
fname = sprintf('0.5 degrees, dyn + interp oldDij.mat');
load(fname)
refDose_YY_oldDij = recalc.resultGUI.physicalDose;

% adjust overlap priorities
cst_Over = matRad_setOverlapPriorities(cst);

V_TargAndNorm = false(size(refDose));
% adjust objectives _and_ constraints internally for fractionation 
for i = 1:size(cst_Over,1)
    for j = 1:size(cst_Over{i,6},1)
       cst_Over{i,6}(j).dose = cst_Over{i,6}(j).dose/pln.numOfFractions;
    end
    if ~isempty(cst_Over{i,6}) && (~strcmp(cst_Over{i,2},'BODY'))
        [x, y, z] = ind2sub(size(refDose),cst_Over{i,4}{1});
        for k = 1:numel(x)
            V_TargAndNorm(x(k),y(k),z(k)) = true;
        end
    end
end

options.lb              = recalc.apertureInfo.limMx(:,1);                                          % Lower bound on the variables.
options.ub              = recalc.apertureInfo.limMx(:,2);                                          % Upper bound on the variables.
[options.cl,options.cu] = matRad_daoGetConstBounds(cst_Over,recalc.apertureInfo,options,recalc.pln.leafSpeedCst,recalc.pln.doseRateCst);   % Lower and upper bounds on the constraint functions.
options.VMAT = recalc.pln.VMAT;

% set optimization options
options.radMod          = recalc.pln.radiationMode;
options.bioOpt          = recalc.pln.bioOptimization;
options.ID              = [recalc.pln.radiationMode '_' recalc.pln.bioOptimization];
%options.numOfScenarios  = dij.numOfScenarios;

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
    dose = recalc.resultGUI.physicalDose;
    percVErr1_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_YY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
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
    dose = recalc.resultGUI.physicalDose;
    percVErr1_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_NY(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_NY(:,:,i) = fluence_NY(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_NY(i) = weight_NY(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NY(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    %next, do interpolation and dynamic fluence, but using the Dij matrices
    %at the original 4degree resolution
    fname = sprintf('%.1f degrees, dyn + interp oldDij.mat',angularRes);
    load(fname);
    dose = recalc.resultGUI.physicalDose;
    percVErr1_YY_oldDij(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_YY_oldDij(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_YY_oldDij(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_YY_oldDij(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    
    percVErr1_YY_oldDij_itself(i) = 100*nnz(abs(dose-refDose_YY_oldDij)./refDose_YY_oldDij >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_YY_oldDij_itself(i) = 100*nnz(abs(dose-refDose_YY_oldDij)./refDose_YY_oldDij >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_YY_oldDij_itself(i) = 100*nnz(abs(dose-refDose_YY_oldDij)./refDose_YY_oldDij >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_YY_oldDij_itself(i) = 100*nnz(abs(dose-refDose_YY_oldDij)./refDose_YY_oldDij >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_YY_oldDij(:,:,i) = fluence_YY_oldDij(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_YY_oldDij(i) = weight_YY_oldDij(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NY(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    %finally, do neither interpolation nor dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + Ninterp.mat',angularRes);
    load(fname);
    dose = recalc.resultGUI.physicalDose;
    percVErr1_NN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_NN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_NN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_NN(i) = 100*nnz(abs(dose-refDose)./refDose >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    
    percVErr1_NN_itself(i) = 100*nnz(abs(dose-refDose_NN)./refDose_NN >= 0.01 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr3_NN_itself(i) = 100*nnz(abs(dose-refDose_NN)./refDose_NN >= 0.03 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr5_NN_itself(i) = 100*nnz(abs(dose-refDose_NN)./refDose_NN >= 0.05 & V_TargAndNorm)./nnz(V_TargAndNorm);
    percVErr10_NN_itself(i) = 100*nnz(abs(dose-refDose_NN)./refDose_NN >= 0.10 & V_TargAndNorm)./nnz(V_TargAndNorm);
    for j = 1:numel(recalc.apertureInfo.beam)
        fluence_NN(:,:,i) = fluence_NN(:,:,i)+recalc.apertureInfo.beam(j).shape(1).weight*recalc.apertureInfo.beam(j).shape(1).shapeMap;
        weight_NN(i) = weight_NN(i)+recalc.apertureInfo.beam(j).shape(1).weight;
    end
    %obj_NN(i) = matRad_daoObjFunc(recalc.apertureInfo.apertureVector,recalc.apertureInfo,dij,cst_Over,options);
    
    
    i = i+1;
end


save('Results','*_NN', '*_NY', '*_YN', '*_YY')

figure
hold
plot(angularResS,percVErr1_YY)
plot(angularResS,percVErr3_YY)
plot(angularResS,percVErr5_YY)
plot(angularResS,percVErr10_YY)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_YY_oldDij)
plot(angularResS,percVErr3_YY_oldDij)
plot(angularResS,percVErr5_YY_oldDij)
plot(angularResS,percVErr10_YY_oldDij)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated, oldDij';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_YY_oldDij_itself)
plot(angularResS,percVErr3_YY_oldDij_itself)
plot(angularResS,percVErr5_YY_oldDij_itself)
plot(angularResS,percVErr10_YY_oldDij_itself)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Dynamic, interpolated, oldDij cf itself';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NY)
plot(angularResS,percVErr3_NY)
plot(angularResS,percVErr5_NY)
plot(angularResS,percVErr10_NY)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NN)
plot(angularResS,percVErr3_NN)
plot(angularResS,percVErr5_NN)
plot(angularResS,percVErr10_NN)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, not interpolated';
title(fname)
grid
savefig(fname)

figure
hold
plot(angularResS,percVErr1_NN_itself)
plot(angularResS,percVErr3_NN_itself)
plot(angularResS,percVErr5_NN_itself)
plot(angularResS,percVErr10_NN_itself)
xlabel('angular resolution (^\circ)')
ylabel('Volume (%)')
legend({'Error > 1%' 'Error > 3%' 'Error > 5%' 'Error > 10%'},'location','best')
fname = 'Not dynamic, not interpolated cf itself';
title(fname)
grid
savefig(fname)



