function [cst] = matRad_genCst(ct,path)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call matRad_genCst
% this function reads out a *.cst file into a cell array for matRad
% 1: organ number in VOI
% 2: organ name
% 3: organ class
% 4: maximum dose
% 5: minimum dose
% 6: maximum penalty
% 7: minimum penalty
% 8: VOILIST for countour.
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
%% Generate CST column 3 to 7

% This columns are handmade generated for every organ.
% You can change this values in PATIENT_NAME_CST_Input files

cst = eval(sprintf(strcat(path,'_CST_Input')));

%% open file for reading
% Save in a struct all files inside folder who finished with "_VOILIST.mat"
files = dir(fullfile(path,'*_VOILIST.mat'));

% Save in a cell all files names.
allNames={files.name};

%% Generate two firsts CST file columns
% This for cicle, runs over every file name and generate the two first
% columns of CST file.
for i = 1:numel(files)
    
    % First column: Organ number. Note: The first number is zero.
    cst{i,1}=i-1;
    
    % Second column: Organ name
    
    % Save names in a temporary variable 'f'.
    f=allNames{i};
    load(f);
    
    % Save names on second column in CST file.
    % Note: The command strrep is used for quit '_VOILIST.mat' characters. Then
    % is used again for replacing '_' character for a space. This is needed
    % for leyend visualization on matRad.
    cst{i,2}=strrep(strrep(f,'_VOILIST.mat',''), '_', ' ');
    
end


%% Generate VOI list. It is storage in column 8.

% This for cicle, runs over every VOILIST organ and save it in a 'V' cell.

for i=1:numel(files)
    load(allNames{i});
    V{i} = v;
end

% Save the VOI list with anatomical (LPS) coordinate system.

% We use the most important model for coordinate system for medical
% imaging, the original data do not use this system, it is only consistent
% for X and Y axis. It is needed to mirror every slice in Z axis.
% More information about LPS system: http://goo.gl/gPM9YP

% This for cicle, runs over every VOILIST reflecting it. For example for a
% Cube who has 90 slices, the slice z=89, is reflected to slice z=2 and
% viceversa.

for s=1:length(V)
    
    % convert Indices to subscripts for every VOILIST
    [y,x,z] = ind2sub(size(ct),V{s});
    
    % Mirror Z position
    z = size(ct,3)-z+1;
    
    % Convert back: subcripts to indices.
    V{s} = sub2ind(size(ct),y,x,z);
    
    % Sort reflected VOILIST and then stored in 8th CST file column.
    cst{s,8} = sort(V{s});
    
end