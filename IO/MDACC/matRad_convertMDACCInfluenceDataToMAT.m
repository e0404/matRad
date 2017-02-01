function [ijFileNames] = matRad_convertMDACCInfluenceDataToMAT(folderName,extension)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

quantitiesToParse = {'*Dij*'}; % {'*Dij*','*LETij*'};
ijCnt             = 1;


for i = 1:numel(quantitiesToParse)
    
files = dir( [folderName filesep quantitiesToParse{1,i} extension]); % get files matching pattern
fileNames = {files.name};                                            % get only file names


    for j = 1:numel(fileNames)

        switch extension
            case {'ASCII'}

                fileName          = [folderName filesep fileNames{1,j}];
                [~,name,~]        = fileparts(fileName);  
                
                [influenceStruct] = readASCII(fileName);
              
                % save influenceStruct to disk
                save([folderName filesep name],'influenceStruct');
                ijFileNames{ijCnt,1} = [folderName filesep name '.mat'];
                clear 'influenceStruct';
                if ~isempty(strfind(name,'Dij'))
                    ijFileNames{ijCnt,2} = 'Dij';
                elseif ~isempty(strfind(name,'LETij'))
                    ijFileNames{ijCnt,2} = 'LETij';
                end
                ijCnt = ijCnt + 1;

            case {'bianry'}
                % implement binary import once corresponding is available
        end

    end
    
end

   
end



function [influceStruct] = readASCII(absolutePath)

blockSize  = 100000;                    % parse in data in blocks
estmSize   = 900e6;                     % maximal number of elements
mSubScript = zeros(estmSize,3,'uint32');   % initialize 
ixBeam     = zeros(estmSize,1,'uint32');
ixBeamlet  = zeros(estmSize,1,'uint32');
data       = zeros(estmSize,1);

fileID     = fopen(absolutePath);

formatSpec = '%.0d%.0d%.0d%.0d%.0d %f';  % format in files is:   119.93.207.0.1000 4.064967e-01

Cnt = 0;

while ~feof(fileID)  
    
    Cnt = Cnt + 1;
	rawData = textscan(fileID,formatSpec,blockSize,'Delimiter','.\t');
    
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
  % TO DO
  absolutePath  = 0; 
  influceStruct = 0;
end