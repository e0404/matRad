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
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0,'Please wait...');
%h.WindowStyle = 'Modal';
steps = 2;

%% import ct-cube
waitbar(1 / steps)
resolution = [files.resx, files.resy, files.resz]; % [mm] / lps coordinate system
ct = matRad_importDicomCt(files.ct, resolution); 
    
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

%% save ct and cst
matRadFileName = [files.rtss{3} '.mat']; % use default from dicom
[FileName,PathName] = uiputfile('*','Save as...',matRadFileName);
if ischar(FileName)
    save([PathName, FileName],'ct','cst');
end
end