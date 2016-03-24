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
ct.nScen = length(InputData);
for i = 1:ct.nScen
    [ct.cube(:,:,:,i), CTcubeReadData3Dinfo(i)] = ReadData3D(fullfile(InputFolder,InputData(i).name,InputData(i).CTs),false);
    display(['import CT ',num2str(i),'/',num2str(ct.nScen)])
end

% HU eDEns conversion
ct.cube = double(ct.cube);
ct.cube = matRad_convHU2eDens(ct.cube);

% set CT resolution
dim = reshape([CTcubeReadData3Dinfo.PixelDimensions],3,size(ct.cube,4))';
if length(unique(dim(:,1))) == 1 && length(unique(dim(:,2))) == 1 && length(unique(dim(:,3))) == 1
    ct.resolution.x = unique(dim(:,1));
    ct.resolution.y = unique(dim(:,2));
    ct.resolution.z = unique(dim(:,3));
else
    error('found CTs with different resolutions');
end

clear CTcubeReadData3Dinfo dim i 

%% read segmentations
display('start segmentation import:')

% check if all VOIs are segmented in all Scenarios
logVOI = true(1,length(VOIs));
for i = 1:ct.nScen
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
    for j = 1:ct.nScen
            idx = find(strcmp(InputData(j).SegmentationsName,VOIs{i}));
            [tmp,~] = ReadData3D(fullfile(InputFolder,InputData(j).name,InputData(j).SegmentationsFile{idx}),false);
            cst{i,4} = [cst{i,4};find(tmp>0)+i*size(ct.cube,1)*size(ct.cube,2)*size(ct.cube,3)];
            display(['import segmentation of VOI ',num2str(i),'/',num2str(length(VOIs)),' in Scenario ',num2str(j),'/',num2str(ct.nScen)])
    end
    
    % Tissue parameters and Dose Objectives
    cst{i,5} = []; 
    cst{i,6} = [];
    
end

display('import finished')
end