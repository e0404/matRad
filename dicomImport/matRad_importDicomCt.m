function ct = matRad_importDicomCt(ctList, resolution, visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to import dicom ct data
% 
% call
%   ct = matRad_importDicomCt(ctList, resolution, visBool)
%
% input
%   ctList:         list of dicom ct files
%   resolution:   	resolution of the imported ct cube, i.e. this function
%                   will interpolate to a different resolution if desired
%   visBool:        optional: turn on/off visualization
%
% output
%   ct:             matRad ct struct. Note that this 3D matlab array 
%                   contains water euqivalen electron denisities.
%                   Hounsfield units are converted using a standard lookup
%                   table in matRad_calcWaterEqD
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

fprintf('\nimporting ct-cube...');

%% processing input variables
if nargin < 3
    visBool = 0;
end

% creation of info list
numOfSlices = size(ctList,1);
fprintf('\ncreating info...')
for i = 1:numOfSlices
    tmpDicomInfo = dicominfo(ctList{i,1});
    
    % remember relevant dicom info - do not record everything as some tags
    % might not been defined for individual files
    info(i).PixelSpacing            = tmpDicomInfo.PixelSpacing;
    info(i).ImagePositionPatient    = tmpDicomInfo.ImagePositionPatient;
    info(i).SliceThickness          = tmpDicomInfo.SliceThickness;
    info(i).ImageOrientationPatient = tmpDicomInfo.ImageOrientationPatient;
    info(i).PatientPosition         = tmpDicomInfo.PatientPosition;
    info(i).Rows                    = tmpDicomInfo.Rows;
    info(i).Columns                 = tmpDicomInfo.Columns;
    info(i).Width                   = tmpDicomInfo.Width;
    info(i).Height                  = tmpDicomInfo.Height;
    info(i).RescaleSlope            = tmpDicomInfo.RescaleSlope;
    info(i).RescaleIntercept        = tmpDicomInfo.RescaleIntercept;
    
    matRad_progress(i,numOfSlices);
end

% adjusting sequence of slices (filenames may not be ordered propperly....
% e.g. CT1.dcm, CT10.dcm, CT100zCoordList = [info.ImagePositionPatient(1,3)]';.dcm, CT101.dcm,...
CoordList = [info.ImagePositionPatient]';
[~, indexing] = sort(CoordList(:,3)); % get sortation from z-coordinates

ctList = ctList(indexing);
info = info(indexing);

%% check data set for consistency
if size(unique([info.PixelSpacing]','rows'),1) > 1
    error('Different pixel size in different CT slices');
end

coordsOfFirstPixel = [info.ImagePositionPatient];
if numel(unique(coordsOfFirstPixel(1,:))) > 1 || numel(unique(coordsOfFirstPixel(2,:))) > 1
    error('Ct slices are not aligned');
end
if sum(diff(coordsOfFirstPixel(3,:))<=0) > 0
    error('Ct slices not monotonically increasing');
end
if numel(unique([info.Rows])) > 1 || numel(unique([info.Columns])) > 1
    error('Ct slice sizes inconsistent');
end


%% checking the patient position
% As of now, the matRad treatment planning system is only valid for
% patients in a supine position. Other orientations (e.g. prone, decubitus
% left/right) are not supported.
% Defined Terms:
% HFP     Head First-Prone                  (not supported)
% HFS     Head First-Supine                 (supported)
% HFDR    Head First-Decubitus Right        (not supported)
% HFDL    Head First-Decubitus Left         (not supported)
% FFDR    Feet First-Decubitus Right        (not supported)
% FFDL    Feet First-Decubitus Left         (not supported)
% FFP     Feet First-Prone                  (not supported)
% FFS     Feet First-Supine                 (supported)

if isempty(regexp(info(1).PatientPosition,'S', 'once'))
    error(['This Patient Position is not supported by matRad.'...
        ' As of now only ''HFS'' (Head First-Supine) and ''FFS'''...
        ' (Feet First-Supine) can be processed.'])    
end

%% creation of ct-cube
fprintf('reading slices...')
origCt = zeros(info(1).Height, info(1).Width, numOfSlices);
for i = 1:numOfSlices
    currentFilename = ctList{i};
    [currentImage, map] = dicomread(currentFilename);
    origCt(:,:,i) = currentImage(:,:); % creation of the ct cube
    
    % draw current ct-slice
    if visBool
        if ~isempty(map)
            image(ind2rgb(uint8(63*currentImage/max(currentImage(:))),map));
            xlabel('x [voxelnumber]')
            ylabel('y [voxelnumber]')
            title(['Slice # ' int2str(i) ' of ' int2str(numOfSlices)])
        else
            image(ind2rgb(uint8(63*currentImage/max(currentImage(:))),bone));
            xlabel('x [voxelnumber]')
            ylabel('y [voxelnumber]')
            title(['Slice # ' int2str(i) ' of ' int2str(numOfSlices)])
        end
        axis equal tight;
        pause(0.1);
    end
    matRad_progress(i,numOfSlices);
end

%% correction if not lps-coordinate-system
% when using the physical coordinates (ctInfo.ImagePositionPatient) to
% arrange the  slices in z-direction, there is no more need for mirroring
% in the z-direction
fprintf('\nz-coordinates taken from ImagePositionPatient\n')

% The x- & y-direction in lps-coordinates are specified in:
% ImageOrientationPatient
xDir = info(1).ImageOrientationPatient(1:3); % lps: [1;0;0]
yDir = info(1).ImageOrientationPatient(4:6); % lps: [0;1;0]
nonStandardDirection = false;

% correct x- & y-direction

if xDir(1) == 1 && xDir(2) == 0 && xDir(3) == 0
    fprintf('x-direction OK\n')
elseif xDir(1) == -1 && xDir(2) == 0 && xDir(3) == 0
    fprintf('\nMirroring x-direction...')
    origCt = flip(origCt,1);
    fprintf('finished!\n')
else
    nonStandardDirection = true;
end
    
if yDir(1) == 0 && yDir(2) == 1 && yDir(3) == 0
    fprintf('y-direction OK\n')
elseif yDir(1) == 0 && yDir(2) == -1 && yDir(3) == 0
    fprintf('\nMirroring y-direction...')
    origCt = flip(origCt,2);
    fprintf('finished!\n')
else
    nonStandardDirection = true;
end

if nonStandardDirection
    fprintf(['Non-standard patient orientation.\n'...
        'CT might not fit to contoured structures\n'])
end

%% interpolate cube
fprintf('\nInterpolating CT cube...');
ct = matRad_interpCtCube(origCt, info, resolution);
fprintf('finished!\n');

%% remember some parameters of original dicom
ct.dicomInfo.PixelSpacing            = info(1).PixelSpacing;
                                       tmp = [info.ImagePositionPatient];
ct.dicomInfo.SlicePositions          = tmp(3,:);
ct.dicomInfo.SliceThickness          = [info.SliceThickness];
ct.dicomInfo.ImagePositionPatient    = info(1).ImagePositionPatient;
ct.dicomInfo.ImageOrientationPatient = info(1).ImageOrientationPatient;
ct.dicomInfo.PatientPosition         = info(1).PatientPosition;
ct.dicomInfo.Width                   = info(1).Width;
ct.dicomInfo.Height                  = info(1).Height;
ct.dicomInfo.RescaleSlope            = info(1).RescaleSlope;
ct.dicomInfo.RescaleIntercept        = info(1).RescaleIntercept;

% convert to water equivalent electron densities
fprintf('\nconversion of ct-Cube to waterEqD...');
ct = matRad_calcWaterEqD(ct);
fprintf('finished!\n');

end