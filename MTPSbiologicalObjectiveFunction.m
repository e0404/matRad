function ObjVal = MTPSbiologicalObjectiveFunction(Weight, DIJ, cst)
% A phyisical objective function

d=DIJ*Weight;
d_i=d(cst{2,8});
d_t=ones(size(d_i))*cst{2,4};


%ObjVal = (AlphaIJ .* DIJ) * Weight + (sqrt(BetaIJ) .* DIJ * Weight).^2  -  PrescribedEffect;
ObjVal =  d_i -  d_t;


%ObjVal(ObjVal(cst{2,8}) < 0) = 0;
ObjVal = (100* ObjVal)' * ObjVal;