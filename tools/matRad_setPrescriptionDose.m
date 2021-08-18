function [cst,pln] = matRad_setPrescriptionDose(prescriptionDose,numOfFractions,cst,pln)

objectives = [cst{:,6}];
if isfield(objectives,'dose')
    cst = matRad_convertOldCstToNewCstObjectives(cst);
end

pln.numOfFractions = numOfFractions;
pln.prescribedDose = prescriptionDose;

for c = 1:size(cst,1)
    if ~isempty(cst{c,6})
        if contains(cst{c,3},{'TARGET'})
                cst{c,6}{1,1}.parameters{1,1} = prescriptionDose;
        else
                cst{c,6}{1,1}.parameters{1,1} = prescriptionDose/12;
        end
    end
end


end