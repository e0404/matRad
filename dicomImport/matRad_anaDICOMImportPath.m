function files = matRad_anaDICOMImportPath(patDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call [ctList, structList, structPath, otherFiles] = analyzePatDir(patDir)
% to analyze the chosen directory ('patDir') and its subdirectory for 
% DICOM files. The output contains a list of all identified CT and
% structure files and the path to a single RTSTRUCT file that shall be
% used for the patient import.
% Unidentified (or files we are not interested in) DICOM files are
% stored in 'otherFiles'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\nAnalyzing patient directory...\n']);

%% get all files in search directory

% get information about main directory
mainDirInfo = dir(patDir);   
% get index of subfolders
dirIndex = [mainDirInfo.isdir]; 
% list of filenames in main directory
fileList = {mainDirInfo(~dirIndex).name}';

% create full path for all files in main directory
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(patDir,x),...
                fileList, 'UniformOutput', false);
end

if isempty(fileList)
    error('Search folder empty')
end

%% check for dicom files and differentiate patients, types, and series
numOfFiles = numel(fileList(:,1));
for i = numOfFiles:-1:1
    try % try to get DicomInfo
        info = dicominfo(fileList{i});
        fileList{i,2} = info.Modality;
        fileList{i,3} = info.PatientID;
        fileList{i,4} = info.SeriesNumber;
    catch
        fileList(i,:) = [];
    end
    matRad_progress(numOfFiles+1-i, numOfFiles);
end

if isempty(fileList)
    error('No DICOM files found in patient directory');
end

%% select patient
patientList = unique(fileList(:,3));

if numel(patientList) > 1
    fprintf('\nFound data sets for multiple patients\n');
    for i = 1:numel(patientList)
        fprintf([ '[' num2str(i) '] ' patientList{i} '\n']);
    end
    patientIx = input(['\nPlease select a patient by entering 1,2,...\n']);
    if isempty(intersect(patientIx,1:numel(patientList)))
        error('Invalid patient identifier\n');
    end
else
    patientIx = 1;
end
fprintf(['Importing patient ' patientList{patientIx} '\n']);

%% select image series
patientFileIx = strcmp(fileList(:,3),patientList{patientIx});
ctFileIx      = strcmp(fileList(:,2),'CT');
rtssFileIx    = strcmp(fileList(:,2),'RTSTRUCT');
       
seriesWithCt   = unique(cell2mat(fileList(ctFileIx&patientFileIx,4)));
seriesWithRtss = unique(cell2mat(fileList(rtssFileIx&patientFileIx,4)));

seriesWithCtAndRtss = intersect(seriesWithCt,seriesWithRtss);

if numel(seriesWithCtAndRtss) < 1
    error('Could not find image series with corresponding rt structure set for selected patient\n');
elseif numel(seriesWithCtAndRtss) == 0
    seriesIx = 1;
    fprintf(['\nImporting image series ' num2str(seriesWithCtAndRtss) ' including corresponding rt structure set\n']);
else
    fprintf('\nFound multiple image series with corresponding rt structure set\n');
    for i = 1:numel(seriesWithCtAndRtss)
        fprintf([ 'Image series [' num2str(seriesWithCtAndRtss(i)) ']\n']);
    end
    seriesIx = input(['\nPlease select an image series by entering the corresponding identifier\n']);
    if isempty(intersect(seriesIx,seriesWithCtAndRtss))
        error('Invalid series identifier\n');
    end
end

%% select rt structure set
patientRtssFileIx =  find(rtssFileIx&patientFileIx&cell2mat(fileList(:,4))==seriesIx);

if isempty(patientRtssFileIx)
    error('Could not find rt structure set for corresponding image series');
else
    if numel(patientRtssFileIx) > 1
        fprintf('\nFound multiple rt structure sets for image series\n');
        for i = 1:numel(patientRtssFileIx)
             display(sprintf('[%d] %s',i,fileList{patientRtssFileIx(i),1}));
        end
        rtssIx = input(['\nPlease select an rt structure set by entering the corresponding identifier\n']);
        if isempty(intersect(rtssIx,1:numel(patientRtssFileIx)))
            error('Invalid series identifier\n');
        end
    end
end

%% remove files that should not be imported from list

files.ct   = fileList(find(ctFileIx&patientFileIx&cell2mat(fileList(:,4))==seriesIx),:);
files.rtss = fileList(rtssIx,:); 


