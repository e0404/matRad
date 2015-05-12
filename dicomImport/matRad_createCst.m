function cst = matRad_createCst(structures)

cst = cell(size(structures,2),6);

for i = 1:size(structures,2)
    cst{i,1} = i - 1; % first organ has number 0    
    cst{i,2} = structures(i).structName;
    
    if ~isempty(regexpi(cst{i,2},'tv')) || ...
       ~isempty(regexpi(cst{i,2},'target')) || ...
       ~isempty(regexpi(cst{i,2},'gtv')) || ...
       ~isempty(regexpi(cst{i,2},'ctv')) || ...
       ~isempty(regexpi(cst{i,2},'ptv')) || ...
       ~isempty(regexpi(cst{i,2},'boost')) || ...
       ~isempty(regexpi(cst{i,2},'tumor'))
        cst{i,3} = 'TARGET';
    else
        cst{i,3} = 'OAR';
    end
    
    cst{i,4} = structures(i).indices;
    cst{i,5}.Priority = 1; % set dummy priority but no biology parameters
    cst{i,6} = []; % define no dummy objcetives    
end
