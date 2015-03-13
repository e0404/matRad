function ObjVal = MTPSbiologicalObjectiveFunction(Weight, DIJ, cst)
% A phyisical objective function
numVoxels = size(DIJ,1);
PrescribedEffect = zeros(numVoxels,1);
PrescribedDose = zeros(numVoxels,1);

PrescribedDose(cst{2,8})=cst{2,4};
PrescribedEffect(cst{2,8}) = 0.1*cst{2,4}+0.05*cst{2,4}^2;

%ObjVal = (AlphaIJ .* DIJ) * Weight + (sqrt(BetaIJ) .* DIJ * Weight).^2  -  PrescribedEffect;
ObjVal =  DIJ*Weight -  PrescribedDose;


%ObjVal(ObjVal(cst{2,8}) < 0) = 0;
ObjVal = (100* ObjVal)' * ObjVal;