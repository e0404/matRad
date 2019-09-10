function [cst] = matRad_slimCST(cst)

j = 1;
while j <= length(cst)
    if contains(cst{j,2},'h_')
        cst(j,:)=[];
        j = j-1;
    elseif contains(cst{j,2},'h.')
        cst(j,:)=[];
        j = j-1;
%     elseif contains(cst{j,2},'CT-')
%         cst(j,:)=[];
    end
    if contains(cst{j,2},'sophagus')
        cst{j,2} = 'Esophagus';
    end
    
    j = j+1;  
end
end