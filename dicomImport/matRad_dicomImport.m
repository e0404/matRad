%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import script to convert DICOM files into matRad base data (i.e. ct,cst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% scan dicom import folder and select patient, image series, and rtss
dicomDataDir = 'X:\YourDicomFolder';
files = matRad_scanDicomImportFolder(dicomDataDir);

%% import ct-cube
resolution = [2 2 2]; % [mm] / lps coordinate system
ct = matRad_importDicomCt(files.ct, resolution); 
    
%% import structure data
structures = matRad_importDicomRtss(files.rtss{1});

%% creating structure cube
for i = 1:numel(structures)
    fprintf('creating cube for %s volume...\n', structures(i).structName);
    structures(i).indices = matRad_convRtssContours2Indices(structures(i).points,ct);
end
fprintf('finished!\n');

%% creating cst
cst = matRad_createCst(structures);

%% save ct and cst
matRadFileName = files.rtss{3}; % use default from dicom
save([matRadFileName '.mat'],'ct','cst');

%% visualize results

% matRad visualization (proper visualization)
pln.resolution = ct.resolution;
pln.isoCenter  = matRad_getIsoCenter(cst,ct);

matRad_visCtDose([],cst,pln,ct);
