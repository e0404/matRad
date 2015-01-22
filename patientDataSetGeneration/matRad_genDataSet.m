function matRad_genDataSet

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call matRad_genMbPlanDataSet
% to to convert a MGH data set into a matlab data set for matRad
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

% Ask for patient to be generated.
result = input('Select patient \n \n 1: Prostate \n 2: Liver \n 3. Head and neck \n 4: TG119 Phantom \n ');

% Set patient path folder and CT resolution [mm/voxel]
if result == 1
    path = 'PROSTATE';
    ctResolution = [3 3 3];
elseif result == 2
    path = 'LIVER';
    ctResolution = [3 3 2.5];
elseif result == 3
    path = 'HN_withoutDij';
    ctResolution = [3 3 5];
elseif result == 4
    path = 'TG119';
    ctResolution = [3 3 2.5];
end

addpath(path);

% Generate the CT file for matRad from MGH data set.
ct = matRad_genCT(result);
disp('CT generated.');

% Generate CST file for matRad from MGH data set.
cst = matRad_genCst(ct,path);
disp('CST file generated.');

% Save CST and CT into .mat file
if result == 1
    save('PROSTATE','cst','ct','ctResolution');
elseif result == 2 
    save('LIVER','cst','ct','ctResolution');
elseif result == 3
    save('HEAD_AND_NECK','cst','ct','ctResolution');
elseif result == 4
    save('TG119','cst','ct','ctResolution');
end

% Note: Strcat is for joint two phrases.
disp(strcat(path,' patient for matRad is already generated.'));

return;