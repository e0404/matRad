%create instance of dose calculation class input: machine basedata. 
% For required fields and their format see matRad/basedata/brachy_LDR or
% class def
isotope = Source(dij.basedata.data); 

switch pln.propDoseCalc.TG43approximation
    case '1D'
        Dose = isotope.getDoseRate1D(DistanceVector);
    case '2D'
        Dose = isotope.getDoseRate2D(DistanceVector,ThetaVector);
end
        
dij.physicalDose = {sparse(reshape(Dose,size(DistanceMatrix.dist)))};
end