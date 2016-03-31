function [ct,cst] = matRad_multScenImport(InputFolder,numOfScen,VOIs)

% add readData3d to search path
addpath('E:\Mescher\12_4DCT\ReadData3d');

%% get input data info
display('start multiple scenario import:')

InputData = dir(InputFolder);
InputData = InputData(~(strcmp('.', {InputData.name}) | strcmp('..', {InputData.name})));
InputData = InputData(1:(min(length(InputData),numOfScen)));

for i = 1:length(InputData)
    InputData(i).CTs = dir(fullfile(InputFolder,InputData(i).name));
    InputData(i).CTs = InputData(i).CTs(~cellfun('isempty',strfind({InputData(i).CTs.name},'_CT_'))); 
    InputData(i).CTs = InputData(i).CTs.name;
    
    InputData(i).SegmentationsFile = dir(fullfile(InputFolder,InputData(i).name,'*.vtk'));
    InputData(i).SegmentationsFile = InputData(i).SegmentationsFile(~cellfun('isempty',strfind({InputData(i).SegmentationsFile.name},'rtss'))); 
    InputData(i).SegmentationsFile = {InputData(i).SegmentationsFile.name};
    for j = 1:length(InputData(i).SegmentationsFile)
        idx = regexp(InputData(i).SegmentationsFile{j},'_');
        InputData(i).SegmentationsName{j} = InputData(i).SegmentationsFile{j}(idx(3)+1:idx(end-2)-2);
    end
end

clear i idx j 

%% read CTs
display('start CT import:')
ct.numOfCtScen = length(InputData);
for i = 1:ct.numOfCtScen
    % read files
    [ct.cube{i}, CTcubeReadData3Dinfo(i)] = ReadData3D(fullfile(InputFolder,InputData(i).name,InputData(i).CTs),false);
    
    % HU eDEns conversion
    ct.cube{i} = double(ct.cube{i});
    ct.cube{i} = matRad_convHU2eDens(ct.cube{i});
    
    % swap x and y (matRad standard)
    ct.cube{i} = permute(ct.cube{i},[2,1,3]);
    
    display(['import CT ',num2str(i),'/',num2str(ct.numOfCtScen)])
end

% set CT resolution
dim = reshape([CTcubeReadData3Dinfo.PixelDimensions],3,ct.numOfCtScen)';
if length(unique(dim(:,1))) == 1 && length(unique(dim(:,2))) == 1 && length(unique(dim(:,3))) == 1
    ct.resolution.x = unique(dim(:,1));
    ct.resolution.y = unique(dim(:,2));
    ct.resolution.z = unique(dim(:,3));
else
    error('found CTs with different resolutions');
end

% set CT cube dimension
cubeDim = reshape([CTcubeReadData3Dinfo.Dimensions],3,ct.numOfCtScen)';
if length(unique(cubeDim(:,1))) == 1 && length(unique(cubeDim(:,2))) == 1 && length(unique(cubeDim(:,3))) == 1
    ct.cubeDim(1) = unique(cubeDim(:,2)); % swap also here x and y
    ct.cubeDim(2) = unique(cubeDim(:,1)); % swap also here x and y
    ct.cubeDim(3) = unique(cubeDim(:,3));
else
    error('found CTs with different resolutions');
end

clear CTcubeReadData3Dinfo dim i cubeDim

%% read segmentations
display('start segmentation import:')

% check if all VOIs are segmented in all Scenarios
logVOI = true(1,length(VOIs));
for i = 1:ct.numOfCtScen
    [logVOItmp,~] = ismember(VOIs,InputData(i).SegmentationsName);
    if sum(logVOItmp) < length(VOIs)
        logVOI = logical(logVOI.*logVOItmp);
        warning(['found no segmentation of ',strjoin({VOIs{~logVOItmp}},' and '),' in Scenario ',InputData(i).name])    
    end
end
if sum(~logVOI) > 0
warning(['do not import segmentation of ',strjoin({VOIs{~logVOI}},' and ')])
end
VOIs = VOIs(logVOI);   

% create cst file
cst = cell(length(VOIs),6);
for i = 1:length(VOIs)
    
    % VOI index
    cst{i,1} = i-1;
    
    % VOI name
    cst{i,2} = VOIs{i};
    
    % VOI type
    if ~isempty(strfind(VOIs{i},'GTV')) | ~isempty(strfind(VOIs{i},'CTV')) | ~isempty(strfind(VOIs{i},'PTV'))
        cst{i,3} = 'TARGET';
    else
        cst{i,3} = 'OAR';
    end
    
    % Voxel indices
    for j = 1:ct.numOfCtScen
            idx     = find(strcmp(InputData(j).SegmentationsName,VOIs{i}));
            [tmp,~] = ReadData3D(fullfile(InputFolder,InputData(j).name,InputData(j).SegmentationsFile{idx}),false);
            
            % swap x and y (matRad standard)
            tmp = permute(tmp,[2,1,3]);
            
            %cst{i,4} = [cst{i,4};find(tmp>0)+(j-1)*ct.cubeDim(1)*ct.cubeDim(2)*ct.cubeDim(3)];
            cst{i,4}{j} = find(tmp>0);
            
            display(['import segmentation of VOI ',num2str(i),'/',num2str(length(VOIs)),' in Scenario ',num2str(j),'/',num2str(ct.numOfCtScen)])
    end
    
    % Tissue parameters
    cst{i,5}.TissueClass = 1;
    cst{i,5}.alphaX      = 0.1; 
    cst{i,5}.betaX       = 0.05;
    cst{i,5}.Priority    = NaN; 
    
    % Dose Objectives
    cst{i,6} = [];
    
end

display('import finished')
end