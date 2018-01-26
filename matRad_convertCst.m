function cst2 = matRad_convertCst(cst)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = size(cst);
k = x(1);

cst2 = cst;

for m=1:k 
    if ~isempty(cst{m,6})
        cst2{m,6} = cell(0);
        for n=1:numel(cst{m,6})
            
            if isequal(cst{m,6}(n).type,'square deviation')
                obj = DoseObjectives.matRad_SquaredDeviation;
                    obj.penalty = cst{m,6}(n).penalty;
                    obj.parameters{1} = cst{m,6}(n).dose;
                cst2{m,6}{n} = obj;
                
            elseif isequal(cst{m,6}(n).type,'square overdosing')
                obj = DoseObjectives.matRad_SquaredOverdosing;
                    obj.penalty = cst{m,6}(n).penalty;
                    obj.parameters{1} = cst{m,6}(n).dose;
                cst2{m,6}{n} = obj;
            
            elseif isequal(cst{m,6}(n).type,'square underdosing')
                obj = DoseObjectives.matRad_SquaredUnderdosing;
                    obj.penalty = cst{m,6}(n).penalty;
                    obj.parameters{1} = cst{m,6}(n).dose;
                cst2{m,6}{n} = obj;
                
            elseif isequal(cst{m,6}(n).type,'MinDVH')
                obj = DoseObjectives.matRad_MinDVH;
                    obj.parameters{1} = cst{m,6}(n).dose;
                    obj.parameters{2} = cst{m,6}(n).volume;
                cst2{m,6}{n} = obj;
                
            elseif isequal(cst{m,6}(n).type,'MaxDVH')
                obj = DoseObjectives.matRad_MaxDVH;
                    obj.parameters{1} = cst{m,6}(n).dose;
                    obj.parameters{2} = cst{m,6}(n).volume;
                cst2{m,6}{n} = obj;
                
            elseif isequal(cst{m,6}(n).type,'MeanDose')
                obj = DoseObjectives.matRad_MeanDose;
                    obj.parameters{1} = cst{m,6}(n).dose;
                cst2{m,6}{n} = obj;
            
            else
                warndlg('ERROR. Can not read file.','Loading Error');
                break;
            end
        end
    end
end
%ans = cst2;

            
                
           



