function [cst] = matRad_applyOARobjective(cst,prescDose)

OAR.type = 'square overdosing';
OAR.penalty = 400;
OAR.dose = prescDose*0.7; % Apply 70% of prescribed dose to OARs
OAR.EUD = NaN;
OAR.volume = NaN;
OAR.robustness = "none";

for i = 1:length(cst)
    if contains(cst{i,3},'OAR') & ~contains(cst{i,2},'ZP') & ~contains(cst{i,2},'CT-') & ~contains(cst{i,2},'Ref')
        cst{i,6} = OAR;
     %   cst{i,5}.Priority = 2;
    end
end
end