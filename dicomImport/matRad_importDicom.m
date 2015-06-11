function [ct, cst] = matRad_importDicom( files )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

h = waitbar(0,'Please wait...');
h.WindowStyle = 'Modal';
steps = 2;
%% import ct-cube
waitbar(1 / steps)
resolution = [2 2 2]; % [mm] / lps coordinate system
ct = matRad_importDicomCt(files.ct, resolution); 
    
%% import structure data
waitbar(2 / steps)
structures = matRad_importDicomRtss(files.rtss{1},ct.dicomInfo);
close(h)

%% creating structure cube
h = waitbar(0,'Please wait...');
h.WindowStyle = 'Modal';
steps = numel(structures);
for i = 1:numel(structures)
    % computations take place here
    waitbar(i / steps)
    fprintf('creating cube for %s volume...\n', structures(i).structName);
    structures(i).indices = matRad_convRtssContours2Indices(structures(i).points,ct);
end
fprintf('finished!\n');
close(h)

%% creating cst
cst = matRad_createCst(structures);

%% save ct and cst
matRadFileName = files.rtss{3}; % use default from dicom
save([matRadFileName '.mat'],'ct','cst');

end