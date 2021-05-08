function obj = analyzeFolder(obj,folder,useOnly)
%ANAYLZEFOLDER Summary of this function goes here
%   Detailed explanation goes here
totFiles = [dir([folder,'/*.dcm']);dir([folder,'/*.RTSS']);dir([folder,'/*.RTPL'])];
notIdentifiedFiles =[];



obj.changedFlag = false(5,1);


for i=1:length(totFiles)
    
    if ~isempty(strfind(totFiles(i).name,'CT')) && ~isempty(strfind(totFiles(i).name,'dcm')) 
        obj.CTfiles{end+1} = [folder,'/',totFiles(i).name];
        obj.changedFlag(1) = true;
    elseif ~isempty(strfind(totFiles(i).name,'RS')) && ~isempty(strfind(totFiles(i).name,'dcm'))
        obj.RSfiles{end+1} = [folder,'/',totFiles(i).name];
        obj.changedFlag(2) = true;
    elseif ~isempty(strfind(totFiles(i).name,'US')) && ~isempty(strfind(totFiles(i).name,'dcm'))
        obj.USfiles{end+1} = [folder,'/',totFiles(i).name];
        obj.changedFlag(3) = true;
    elseif ~isempty(strfind(totFiles(i).name,'MR')) && ~isempty(strfind(totFiles(i).name,'dcm'))
        obj.MRfiles{end+1} = [folder,'/',totFiles(i).name];
        obj.changedFlag(4) = true;
    elseif ~isempty(strfind(totFiles(i).name,'RP')) && ~isempty(strfind(totFiles(i).name,'dcm'))
        obj.RPfiles{end+1} = [folder,'/',totFiles(i).name];
        obj.changedFlag(5) = true;
    else
        notIdentifiedFiles{end+1} = totFiles(i).name;
    end
    
end

if ~isempty(notIdentifiedFiles)
    for i=1:length(notIdentifiedFiles)
        
        inf = dicominfo([folder,'/',notIdentifiedFiles{i}]);
        if strcmp(inf.Modality,'CT')
            obj.CTfiles{end+1} = [folder,'/',notIdentifiedFiles{i}];
            obj.changedFlag(1) = true;
        elseif strcmp(inf.Modality,'RTSTRUCT')
            obj.RSfiles{end+1} = [folder,'/',notIdentifiedFiles{i}];
            obj.changedFlag(2) = true;
        elseif strcmp(inf.Modality,'US')
            obj.USfiles{end+1} = [folder,'/',notIdentifiedFiles{i}];
            obj.changedFlag(3) = true;
        elseif strcmp(inf.Modality,'RTPLAN')
            obj.RPfiles{end+1} = [folder,'/',notIdentifiedFiles{i}];
            obj.changedFlag(5) = true;
        elseif strcmp(inf.Modality,'MR')
            obj.MRfiles{end+1} = [folder,'/',notIdentifiedFiles{i}];
            obj.changedFlag(4) = true;
        else
            
            fprintf('Warning: Modality "%s" not supported\n',inf.Modality)
        end
    end
end





end

