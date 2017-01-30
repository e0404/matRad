function [ dijFileNames,LETFileNames ] = matRad_convertMDACCInfluenceDataToMAT(folderName,extension)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% parse firstly dij files
files = dir( [folderName filesep '*Dij*.' extension]); %get files matching pattern
fileNames = { files.name}; %get only file names


for i = 1:numel(fileNames)
    
    switch extension
        case {'ASCII'}
            
            fileName        =[folderName filesep fileNames{1,i}];
            
            [influenceStruct] = readASCII(fileName);
            
            [~,name,~] = fileparts(fileName);       
            save([folderName filesep name],'influenceStruct');
            dijFileNames{i} = [folderName filesep name '.mat'];
            
        case {'bianry'}
            
    end
          
end


%% parse LET files
files = dir( [folderName filesep '*LET*.' extension]); %get files matching pattern
fileNames = { files.name}; %get only file names

for i = 1:numel(fileNames)
    
    switch extension
        case {'ASCII'}
            
            fileName        =[folderName filesep fileNames{1,i}];
            [influenceStruct] = readASCII(fileName);
            [~,name,~] = fileparts(fileName);       
            save([folderName filesep name],'influenceStruct');
            LETFileNames{i} = [folderName filesep name '.mat'];
            
        case {'bianry'}
            
    end
          
end
   
   
end



function [influceStruct] = readASCII(absolutePath)

BlockSize  = 100000;

mSubScript = zeros(500e6,3,'uint32');
ixBeam     = zeros(500e6,1,'uint32');
ixBeamlet  = zeros(500e6,1,'uint32');
data       = zeros(500e6,1);

fileID     = fopen(absolutePath);

formatSpec = '%.0d%.0d%.0d%.0d%.0d %f';  % 119.93.207.0.1000 4.064967e-01

Cnt = 0;

while ~feof(fileID)  
    
    Cnt = Cnt + 1;
	rawData = textscan(fileID,formatSpec,BlockSize,'Delimiter','.\t');
    
    currentLenght  = size(rawData{1,1},1);
    mSubScript(Cnt:Cnt+currentLenght-1,1) = rawData{1,1};
    mSubScript(Cnt:Cnt+currentLenght-1,2) = rawData{1,2};
    mSubScript(Cnt:Cnt+currentLenght-1,3) = rawData{1,3};
    ixBeam(Cnt:Cnt+currentLenght-1,1)     = rawData{1,4};
    ixBeamlet(Cnt:Cnt+currentLenght-1,1)  = rawData{1,5};
    data(Cnt:Cnt+currentLenght-1,1)       = rawData{1,6};
   
    Cnt = Cnt + currentLenght-1; 
    
end

mSubScript = mSubScript(1:Cnt,:);
ixBeam     = ixBeam(1:Cnt,:);
ixBeamlet  = ixBeamlet(1:Cnt,:);
data       = data(1:Cnt,:);

influceStruct.mSubScript       = mSubScript;
influceStruct.stats.maxX       = max(mSubScript(:,1));
influceStruct.stats.maxY       = max(mSubScript(:,2));
influceStruct.stats.maxZ       = max(mSubScript(:,3));
influceStruct.stats.minX       = min(mSubScript(:,1));
influceStruct.stats.minY       = min(mSubScript(:,2));
influceStruct.stats.minZ       = min(mSubScript(:,3));
influceStruct.ixBeam           = ixBeam;
influceStruct.ixBeamlet        = ixBeamlet;
influceStruct.data             = data;


end


function [influceStruct] = readBinary(absolutePath)

end