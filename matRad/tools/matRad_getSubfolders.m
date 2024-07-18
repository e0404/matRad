function [subFolder] = matRad_getSubfolders(folderName)
    
    subFolder = {folderName};
    currSubFolders = dir(folderName);

    % Filter normal files
    currSubFolders = currSubFolders([currSubFolders.isdir]);
    currSubFolders = currSubFolders(~ismember({currSubFolders(:).name}, {'.', '..'}));
        
    for subIdx=1:numel(currSubFolders)
        subFolder = [subFolder, matRad_getSubfolders(fullfile(currSubFolders(subIdx).folder, currSubFolders(subIdx).name))];
    end

end