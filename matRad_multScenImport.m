function [ct,cst] = matRad_multScenImport(InputFolder,numOfScen,VOIs)
InputFolder = 'E:\Mescher\13_BIOM_model\01_BioMechModel\Input';
numOfScen = 3;
VOIs = {'blase','darm','haut','prostata','rektum','restblase','rueckenmark',...
        'CTV','PTV','GTV'};

% add readData3d to search path
addpath('E:\Mescher\12_4DCT\ReadData3d');

%% get input data info
InputData = dir(InputFolder);
InputData = InputData(~(strcmp('.', {InputData.name}) | strcmp('..', {InputData.name})));
InputData = InputData(1:(min(length(InputData),numOfScen)));

for i = 1:length(InputData)
    InputData(i).CTs = dir(fullfile(InputFolder,InputData(i).name));
    InputData(i).CTs = InputData(i).CTs(~cellfun('isempty',strfind({InputData(i).CTs.name},'_CT_'))); 
    InputData(i).CTs = InputData(i).CTs.name;
    
    InputData(i).Segmentations = dir(fullfile(InputFolder,InputData(i).name,'*.vtk'));
    InputData(i).Segmentations = InputData(i).Segmentations(~cellfun('isempty',strfind({InputData(i).Segmentations.name},'rtss'))); 
    InputData(i).Segmentations = {InputData(i).Segmentations.name};
end

%% read CTs
ct.nPhases = length(InputData);
for i = 1:ct.nPhases
    [ct.cube(:,:,:,i), CTcubeReadData3Dinfo(i)] = ReadData3D(fullfile(InputFolder,InputData(i).name,InputData(i).CTs),false);
    display(['import CT ',num2str(i),'/',num2str(ct.nPhases)])
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

%% read segmentations
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
    for j = 1:ct.nPhases
        if ~isempty(strfind(InputData(j).Segmentations,VOIs{i}))
            [tmp,~] = ReadData3D(fullfile(InputFolder,InputData(j).name,InputData.Segmentations),false);
            cst{i,4} = [cst{i,4},find(tmp>0)+i*prod(ct.dimensions)]
        else
            error(['found no segmentation of "',VOIs{i},'" in Scenario ',InputData(j).name])
        end
    end
    
    % Tissue parameters and Dose Obkectives
    cst{i,5} = []; 
    cst{i,6} = [];
    

%     for j = 1:length(segmentations)
%         cst(j,:,i) = {j-1,segmentations(j).name(16:(end-10)),'',[],[],[]};
%         [tmp,~] = ReadData3D(fullfile(patientFolder,'biomech_samples',num2str(i+1),segmentations(j).name), false);
%         cst{j,4,i} = find(tmp>0);
%     end
%     display(['import ',num2str(i),'/',num2str(ndefPhases)])
end


end