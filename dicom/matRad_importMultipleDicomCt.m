function [ct,cst] = matRad_importMultipleDicomCt(files_dir,metadata)
% matRad wrapper function to import and combine a set of dicom files 
% into matRad's native data formats
% 
% call
% [ct, cst] = matRad_importMultipleDicomCt(files)  
% [ct, cst] = matRad_importMultipleDicomCt(files,metadata)
%
% input
%   files:          list of files to be imported (will contain ct and rtss
%                   structure set)
%   metadata:       struct of metadata
%
% output
%   ct:        matRad ct multi-scenario struct
%   cst:       matRad cst multi-scenario  struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[env, ~] = matRad_getEnvironment();

%%
matRad_cfg = MatRad_Config.instance();

% Scanning the container directory
[fileList, patient_listbox] = matRad_scanDicomImportFolder([matRad_cfg.matRadRoot files_dir]);

% Getting info about the CT scenarios 
ct.numOfCtScen = length(patient_listbox);
ct.timeStamp = datestr(clock);
ct.resolution.x = metadata.resolution(1);
ct.resolution.y = metadata.resolution(2);
ct.resolution.z = metadata.resolution(3);

files.resx = ct.resolution.x;
files.resy = ct.resolution.y;
files.resz = ct.resolution.z;

if ~isfield(metadata,'useDoseGrid')
	files.useDoseGrid=false;
end

% Initializing ct and cst struct for each scenario 
tmp_ct = cell(1,ct.numOfCtScen);
tmp_cst = cell(1,ct.numOfCtScen);


for ctPhase = 1:ct.numOfCtScen

    files.ct = unique(fileList(strcmp(fileList(:,2), 'CT') & strcmp(fileList(:,3), patient_listbox{ctPhase})));
    files.rtss = unique(fileList(strcmp(fileList(:,2), 'RTSTRUCT') & strcmp(fileList(:,3), patient_listbox{ctPhase})));
    % Importing the Dicom info into cst and cst struct
    [tmp_ct{ctPhase},tmp_cst{ctPhase},~,~,~]=matRad_importDicom(files, true);
    
end

% Checking dimentions and structures for all scenarios
if(all(cellfun(@(ct) isequal(ct.x,tmp_ct{1}.x), tmp_ct)))
    ct.x=tmp_ct{1}.x;
else
    error('There are scenarios with incompatible value of x-coordinates \n');
end

if(all(cellfun(@(ct) isequal(ct.y,tmp_ct{1}.y), tmp_ct)))
    ct.y=tmp_ct{1}.y;
else
    error('There are scenarios with incompatible value of y-coordinates \n');
end

if(all(cellfun(@(ct) isequal(ct.z,tmp_ct{1}.z), tmp_ct)))
    ct.z=tmp_ct{1}.z;
else
   error('There are scenarios with incompatible value of z-coordinates \n');
end

if(all(cellfun(@(ct) isequal(ct.cubeDim,tmp_ct{1}.cubeDim), tmp_ct)))
    ct.cubeDim=tmp_ct{1}.cubeDim;
else
    error('There are scenarios with incompatible value of cubeDim \n');
end

if(all(cellfun(@(cst) isequal([cst{:,1:3}],[tmp_cst{1}{:,1:3}]), tmp_cst))...
        && all(cellfun(@(cst) isequal([cst{:,5}],[tmp_cst{1}{:,5}]), tmp_cst))...
        && all(cellfun(@(cst) isequal([cst{:,6}],[tmp_cst{1}{:,6}]), tmp_cst)))
    cst=tmp_cst{1};
    [numOfStruct, ~] = size(cst);
else
    error('There are scenarios with incompatible structures data \n');
end


% Combining ct and cst info into a multiscenario struct
for ctPhase = 1:ct.numOfCtScen
    
    ct.cubeHU{ctPhase} = tmp_ct{ctPhase}.cubeHU{1};
    ct.dicomInfo(ctPhase) = tmp_ct{ctPhase}.dicomInfo;
    ct.dicomMeta(ctPhase) = tmp_ct{ctPhase}.dicomMeta;

    for structure = 1:numOfStruct
        cst{structure,4}{1,ctPhase}=tmp_cst{ctPhase}{structure,4}{1}; 
    end
    
end

end

