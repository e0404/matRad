function ct = matRad_importDicomCt(ctList, resolution, dicomMetaBool, grid, visBool)
% matRad function to import dicom ct data
% 
% call
%   ct = matRad_importDicomCt(ctList, resolution, dicomMetaBool)
%   ct = matRad_importDicomCt(ctList, resolution, dicomMetaBool, grid)
%   ct = matRad_importDicomCt(ctList, resolution, dicomMetaBool, visBool)
%   ct = matRad_importDicomCt(ctList, resolution, dicomMetaBool, grid, visBool)
%
% input
%   ctList:         list of dicom ct files
%   resolution:   	resolution of the imported ct cube, i.e. this function
%                   will interpolate to a different resolution if desired
%   dicomMetaBool:  store complete dicom information if true
%   grid:           optional: a priori grid specified for interpolation
%   visBool:        optional: turn on/off visualization
%
% output
%   ct:             matRad ct struct. Note that this 3D matlab array 
%                   contains water euqivalent electron denisities.
%                   Hounsfield units are converted using a standard lookup
%                   table in matRad_calcWaterEqD
%
% References
%   -
%
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

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('\nimporting ct-cube...');

%% processing input variables
if ~exist('visBool','var')
  visBool = 0;
end

% creation of ctInfo list
numOfSlices = size(ctList,1);
matRad_cfg.dispInfo('\ncreating info...')

sliceThicknessStandard = true;
for i = 1:numOfSlices

    if verLessThan('matlab','9')
        tmpDicomInfo = dicominfo(ctList{i,1});
    else
        tmpDicomInfo = dicominfo(ctList{i,1},'UseDictionaryVR',true);
    end
    
    % remember relevant dicom info - do not record everything as some tags
    % might not been defined for individual files
    ctInfo(i).PixelSpacing            = tmpDicomInfo.PixelSpacing;
    ctInfo(i).ImagePositionPatient    = tmpDicomInfo.ImagePositionPatient;
    ctInfo(i).SliceThickness          = tmpDicomInfo.SliceThickness;
    ctInfo(i).ImageOrientationPatient = tmpDicomInfo.ImageOrientationPatient;
    ctInfo(i).PatientPosition         = tmpDicomInfo.PatientPosition;
    ctInfo(i).Rows                    = tmpDicomInfo.Rows;
    ctInfo(i).Columns                 = tmpDicomInfo.Columns;
    ctInfo(i).Width                   = tmpDicomInfo.Width;
    ctInfo(i).Height                  = tmpDicomInfo.Height;
    ctInfo(i).RescaleSlope            = tmpDicomInfo.RescaleSlope;
    ctInfo(i).RescaleIntercept        = tmpDicomInfo.RescaleIntercept;
    
    %Problem due to some CT files using non-standard SpacingBetweenSlices
    
    if isempty(ctInfo(i).SliceThickness)
        %Print warning ocne
        if sliceThicknessStandard
            matRad_cfg.dispWarning('Non-standard use of SliceThickness Attribute (empty), trying to overwrite with SpacingBetweenSlices');
            sliceThicknessStandard = false;
        end
        ctInfo(i).SliceThickness = tmpDicomInfo.SpacingBetweenSlices;
    end
    
    if i == 1
        completeDicom = tmpDicomInfo;
    end
    
    matRad_progress(i,numOfSlices);
end

% adjusting sequence of slices (filenames may not be ordered propperly....
% e.g. CT1.dcm, CT10.dcm, CT100zCoordList = [ctInfo.ImagePositionPatient(1,3)]';.dcm, CT101.dcm,...
CoordList = [ctInfo.ImagePositionPatient]';
[~, indexing] = sort(CoordList(:,3)); % get sortation from z-coordinates

ctList = ctList(indexing);
ctInfo = ctInfo(indexing);

%% check data set for consistency
if size(unique([ctInfo.PixelSpacing]','rows'),1) > 1
    matRad_cfg.dispError('Different pixel size in different CT slices');
end

coordsOfFirstPixel = [ctInfo.ImagePositionPatient];
if numel(unique(coordsOfFirstPixel(1,:))) > 1 || numel(unique(coordsOfFirstPixel(2,:))) > 1
    matRad_cfg.dispError('Ct slices are not aligned');
end
if sum(diff(coordsOfFirstPixel(3,:))<=0) > 0
    matRad_cfg.dispError('Ct slices not monotonically increasing');
end
if numel(unique([ctInfo.Rows])) > 1 || numel(unique([ctInfo.Columns])) > 1
    matRad_cfg.dispError('Ct slice sizes inconsistent');
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

if isempty(regexp(ctInfo(1).PatientPosition,{'S','P'}, 'once'))
    matRad_cfg.dispError(['This Patient Position is not supported by matRad.'...
        ' As of now only ''HFS'' (Head First-Supine), ''FFS'''...
        ' (Feet First-Supine), '...    
        '''HFP'' (Head First-Prone), and ''FFP'''...
        ' (Feet First-Prone) can be processed.'])    
end

%% creation of ct-cube
matRad_cfg.dispInfo('reading slices...')
origCt = zeros(ctInfo(1).Height, ctInfo(1).Width, numOfSlices);
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
matRad_cfg.dispInfo('\nz-coordinates taken from ImagePositionPatient\n')

% The x- & y-direction in lps-coordinates are specified in:
% ImageOrientationPatient
xDir = ctInfo(1).ImageOrientationPatient(1:3); % lps: [1;0;0]
yDir = ctInfo(1).ImageOrientationPatient(4:6); % lps: [0;1;0]
nonStandardDirection = false;

% correct x- & y-direction
% 
% if xDir(1) == 1 && xDir(2) == 0 && xDir(3) == 0
%     matRad_cfg.dispInfo('x-direction OK\n')
% elseif xDir(1) == -1 && xDir(2) == 0 && xDir(3) == 0
%     matRad_cfg.dispInfo('\nMirroring x-direction...')
%     origCt = flip(origCt,1);
%     matRad_cfg.dispInfo('finished!\n')
% else
%     nonStandardDirection = true;
% end
%     
% if yDir(1) == 0 && yDir(2) == 1 && yDir(3) == 0
%     matRad_cfg.dispInfo('y-direction OK\n')
% elseif yDir(1) == 0 && yDir(2) == -1 && yDir(3) == 0
%     matRad_cfg.dispInfo('\nMirroring y-direction...')
%     origCt = flip(origCt,2);
%     matRad_cfg.dispInfo('finished!\n')
% else
%     nonStandardDirection = true;
% end

if nonStandardDirection
    matRad_cfg.dispInfo(['Non-standard patient orientation.\n'...
        'CT might not fit to contoured structures\n'])
end

%% interpolate cube
matRad_cfg.dispInfo('\nInterpolating CT cube...');
if exist('grid','var')
    ct = matRad_interpDicomCtCube(origCt, ctInfo, resolution, grid);
else
    ct = matRad_interpDicomCtCube(origCt, ctInfo, resolution);
end
matRad_cfg.dispInfo('finished!\n');

%% remember some parameters of original dicom
ct.dicomInfo.PixelSpacing            = ctInfo(1).PixelSpacing;
                                       tmp = [ctInfo.ImagePositionPatient];
ct.dicomInfo.SlicePositions          = tmp(3,:);
ct.dicomInfo.SliceThickness          = [ctInfo.SliceThickness];
ct.dicomInfo.ImagePositionPatient    = ctInfo(1).ImagePositionPatient;
ct.dicomInfo.ImageOrientationPatient = ctInfo(1).ImageOrientationPatient;
ct.dicomInfo.PatientPosition         = ctInfo(1).PatientPosition;
ct.dicomInfo.Width                   = ctInfo(1).Width;
ct.dicomInfo.Height                  = ctInfo(1).Height;
ct.dicomInfo.RescaleSlope            = ctInfo(1).RescaleSlope;
ct.dicomInfo.RescaleIntercept        = ctInfo(1).RescaleIntercept;
if isfield(completeDicom, 'Manufacturer')
ct.dicomInfo.Manufacturer            = completeDicom.Manufacturer;
end
if isfield(completeDicom, 'ManufacturerModelName')
ct.dicomInfo.ManufacturerModelName   = completeDicom.ManufacturerModelName;
end
if isfield(completeDicom, 'ConvolutionKernel')
ct.dicomInfo.ConvolutionKernel       = completeDicom.ConvolutionKernel;
end

% store patientName only if user wants to
if isfield(completeDicom,'PatientName') && dicomMetaBool == true
    ct.dicomInfo.PatientName         = completeDicom.PatientName;
end
if dicomMetaBool == true
    ct.dicomMeta                     = completeDicom;
end

ct.timeStamp = datestr(clock);

% convert to Hounsfield units
matRad_cfg.dispInfo('\nconversion of ct-Cube to Hounsfield units...');
ct = matRad_calcHU(ct);
matRad_cfg.dispInfo('finished!\n');

end
