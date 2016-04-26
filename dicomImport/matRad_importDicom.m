function [ct, cst] = matRad_importDicom( files )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad wrapper function to import a predefined set of dicom files into
% matRad's native data formats
% 
% call
%   [ct, cst] = matRad_importDicom( files )
%
% input
%   files:  list of files to be imported (will contain cts and rt structure set)
%
% output
%   ct:     matRad ct struct
%   cst:    matRad cst struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0,'Please wait...');
%h.WindowStyle = 'Modal';
steps = 2;

%% import ct-cube
waitbar(1 / steps)
resolution.x = files.resx;
resolution.y = files.resy;
resolution.z = files.resz; % [mm] / lps coordinate system
ct = matRad_importDicomCt(files.ct, resolution); 
    
if ~isempty(files.rtss)
    
    %% import structure data
    waitbar(2 / steps)
    structures = matRad_importDicomRtss(files.rtss{1},ct.dicomInfo);
    close(h)

    %% creating structure cube
    h = waitbar(0,'Please wait...');
    %h.WindowStyle = 'Modal';
    steps = numel(structures);
    for i = 1:numel(structures)
        % computations take place here
        waitbar(i / steps)
        fprintf('creating cube for %s volume...\n', structures(i).structName);
        structures(i).indices = matRad_convRtssContours2Indices(structures(i),ct);
    end
    fprintf('finished!\n');
    close(h)

    %% creating cst
    cst = matRad_createCst(structures);

else
    
    cst = matRad_dummyCst(ct);
    
end
%% save ct and cst
matRadFileName = [files.ct{1,3} '.mat']; % use default from dicom
[FileName,PathName] = uiputfile('*','Save as...',matRadFileName);
if ischar(FileName)
    save([PathName, FileName],'ct','cst');
end
end