function [ sData ] = matRadParseBioData(PathToBaseData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Files = dir(PathToBaseData);
Headers = cell(0);
Cnt = 1;
for i=1:length(Files)
    if Files(i).isdir == 0 && Files(i).bytes >0
        
        str = ['C_E' num2str(sscanf(Files(i).name,'C_E%d')) '_'];
        if isempty(find(ismember(Headers,str)))
            FoundFiles = dir([PathToBaseData filesep str '*']);
            Headers{Cnt}=str;
        else
            FoundFiles = [];
        end
        
        for j=1:length(FoundFiles)
            
            if  ~isempty(strfind(FoundFiles(j).name, 'alpha')) && ~isempty(strfind(FoundFiles(j).name,  Headers{Cnt}))
                Import =importdata([PathToBaseData  filesep  FoundFiles(j).name]);
                sData{Cnt}.energy = sscanf(char(Import.textdata(4,1)),'!Energy: %f');
                sData{Cnt}.depths = Import.data(:,1);
                sData{Cnt}.dEdxA = Import.data(:,2);

            end

            if ~isempty(strfind(FoundFiles(j).name, 'beta')) && ~isempty(strfind(FoundFiles(j).name,  Headers{Cnt}))
                 Import =importdata([PathToBaseData  filesep  FoundFiles(j).name]);
                 sData{Cnt}.dEdxB = Import.data(:,2);
            end
            
            if j==length(FoundFiles)
                Cnt=Cnt+1;
            end
        end
    end
end


A= (cell2mat(sData));
vEnergies=[A(:).energy];
[val idx] = sort(vEnergies);

sData= (cell2mat(sData(idx)));

%  figure,hold on,xlabel('depth in cm'), ylabel('Gy-1')
% for i = 1:5:length(sData)
%      str =['dE/dx * alpha with energy: '  num2str(sData(i).energy)];
%     
%     plot(sData(i).depths,sData(i).dEdxA), title(str,'FontSize',18),grid on
%     waitforbuttonpress
% end

% 
%  figure,hold on,xlabel('depth in cm'), ylabel('Gy-1')
% for i = 1:5:length(sData)
%      str =['dE/dx * sqrt(beta) with energy: '  num2str(sData(i).Energy)];
%     
%     plot(sData(i).Depth,sData(i).dEdxB), title(str,'FontSize',18),grid on
%     waitforbuttonpress
% end

end

