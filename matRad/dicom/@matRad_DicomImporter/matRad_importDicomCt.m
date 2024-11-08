function obj = matRad_importDicomCt(obj)
% matRad function to import dicom ct data
% 
% In your object, there must be properties that contain:
%   - list of dicom ct files;
%   - resolution of the imported ct cube, i.e. this function will 
%   interpolate to a different resolution if desired;
%   - a boolean, if you don't want to import complete dicom information set
%   it false.
% Optional:
%   - a priori grid specified for interpolation;
%   - a boolean to turn off/on visualization.
%
% Output - matRad ct structure. 
% Note that this 3D matlab array contains water euqivalent 
% electron denisities. Hounsfield units are converted using a standard
% lookup table in matRad_calcWaterEqD
%
% call
%   matRad_importDicomCt(obj)
%
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('\nimporting ct-cube...');

%% processing input variables
if ~exist('visBool','var')
  obj.visBool = 0;
end

matRad_checkEnvDicomRequirements(matRad_cfg.env);


% creation of ctInfo list
numOfSlices = size(obj.importFiles.ct, 1);
matRad_cfg.dispInfo('\ncreating info...')

sliceThicknessStandard = true;
for i = 1:numOfSlices

    if matRad_cfg.isOctave || verLessThan('matlab','9')
        tmpDicomInfo = dicominfo(obj.importFiles.ct{i, 1});
    else
        tmpDicomInfo = dicominfo(obj.importFiles.ct{i, 1},'UseDictionaryVR',true);
    end
    
    % remember relevant dicom info - do not record everything as some tags
    % might not been defined for individual files
    obj.importCT.ctInfo(i).PixelSpacing            = tmpDicomInfo.PixelSpacing;
    obj.importCT.ctInfo(i).ImagePositionPatient    = tmpDicomInfo.ImagePositionPatient;
    obj.importCT.ctInfo(i).SliceThickness          = str2double(obj.importFiles.resz);
    obj.importCT.ctInfo(i).ImageOrientationPatient = tmpDicomInfo.ImageOrientationPatient;
    obj.importCT.ctInfo(i).PatientPosition         = tmpDicomInfo.PatientPosition;
    obj.importCT.ctInfo(i).Rows                    = tmpDicomInfo.Rows;
    obj.importCT.ctInfo(i).Columns                 = tmpDicomInfo.Columns;
    obj.importCT.ctInfo(i).Width                   = tmpDicomInfo.Columns;%tmpDicomInfo.Width;
    obj.importCT.ctInfo(i).Height                  = tmpDicomInfo.Rows;%tmpDicomInfo.Height;
    obj.importCT.ctInfo(i).RescaleSlope            = tmpDicomInfo.RescaleSlope;
    obj.importCT.ctInfo(i).RescaleIntercept        = tmpDicomInfo.RescaleIntercept;
    
    %Problem due to some CT files using non-standard SpacingBetweenSlices
    
    if isempty(obj.importCT.ctInfo(i).SliceThickness)
        %Print warning once
        if sliceThicknessStandard
            matRad_cfg.dispWarning('Non-standard use of SliceThickness Attribute (empty), trying to overwrite with SpacingBetweenSlices');
            sliceThicknessStandard = false;
        end
        obj.importCT.ctInfo(i).SliceThickness = tmpDicomInfo.SpacingBetweenSlices;
    end
    
    if i == 1
        completeDicom = tmpDicomInfo;
    end
    
    % Show progress
    if matRad_cfg.logLevel > 2
        matRad_progress(i,numOfSlices);
    end
end

% adjusting sequence of slices (filenames may not be ordered propperly....
% e.g. CT1.dcm, CT10.dcm, CT100zCoordList = [ctInfo.ImagePositionPatient(1,3)]';.dcm, CT101.dcm,...
CoordList = [obj.importCT.ctInfo.ImagePositionPatient]';
[~, indexing] = unique(CoordList(:,3)); % get sortation from z-coordinates

obj.importFiles.ct = obj.importFiles.ct(indexing');
obj.importCT.ctInfo = obj.importCT.ctInfo(indexing');

%% check data set for consistency
if size(unique([obj.importCT.ctInfo.PixelSpacing]','rows'),1) > 1
    matRad_cfg.dispError('Different pixel size in different CT slices');
end

coordsOfFirstPixel = [obj.importCT.ctInfo.ImagePositionPatient];
if numel(unique(coordsOfFirstPixel(1,:))) > 1 || numel(unique(coordsOfFirstPixel(2,:))) > 1
    matRad_cfg.dispError('Ct slices are not aligned');
end
if sum(diff(coordsOfFirstPixel(3,:))<=0) > 0
    matRad_cfg.dispError('Ct slices not monotonically increasing');
end
if numel(unique([obj.importCT.ctInfo.Rows])) > 1 || numel(unique([obj.importCT.ctInfo.Columns])) > 1
    matRad_cfg.dispError('Ct slice sizes inconsistent');
end


%% checking the patient position
% As of now, the matRad treatment planning system is only valid for
% patients in a supine position. Other orientations (e.g. prone, decubitus
% left/right) are not supported.
% Defined Terms:
% HFP     Head First-Prone                  (supported)
% HFS     Head First-Supine                 (supported)
% HFDR    Head First-Decubitus Right        (not supported)
% HFDL    Head First-Decubitus Left         (not supported)
% FFDR    Feet First-Decubitus Right        (not supported)
% FFDL    Feet First-Decubitus Left         (not supported)
% FFP     Feet First-Prone                  (supported)
% FFS     Feet First-Supine                 (supported)

if isempty(regexp(obj.importCT.ctInfo(1).PatientPosition,{'S','P'}, 'once'))
    matRad_cfg.dispError(['This Patient Position is not supported by matRad.'...
        ' As of now only ''HFS'' (Head First-Supine), ''FFS'''...
        ' (Feet First-Supine), '...    
        '''HFP'' (Head First-Prone), and ''FFP'''...
        ' (Feet First-Prone) can be processed.'])    
end

%% creation of ct-cube
matRad_cfg.dispInfo('reading slices...')
origCt = zeros(obj.importCT.ctInfo(1).Height, obj.importCT.ctInfo(1).Width, numOfSlices);
for i = 1:numOfSlices
    currentFilename = obj.importFiles.ct{i};
    if matRad_cfg.isOctave
        currentImage = dicomread(currentFilename);
        map = [];
    else
        [currentImage, map] = dicomread(currentFilename);
    end    
    origCt(:,:,i) = currentImage(:,:); % creation of the ct cube
    
    % draw current ct-slice
    if obj.visBool
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
    % Show progress
    if matRad_cfg.logLevel > 2
        matRad_progress(i,numOfSlices);
    end
end

%% correction if not lps-coordinate-system
% when using the physical coordinates (ctInfo.ImagePositionPatient) to
% arrange the  slices in z-direction, there is no more need for mirroring
% in the z-direction
matRad_cfg.dispInfo('\nz-coordinates taken from ImagePositionPatient\n')

% The x- & y-direction in lps-coordinates are specified in:
% ImageOrientationPatient
xDir = obj.importCT.ctInfo(1).ImageOrientationPatient(1:3); % lps: [1;0;0]
yDir = obj.importCT.ctInfo(1).ImageOrientationPatient(4:6); % lps: [0;1;0]
nonStandardDirection = false;

% correct x- & y-direction

if xDir(1) == 1 && xDir(2) == 0 && xDir(3) == 0
    matRad_cfg.dispInfo('x-direction OK\n')
elseif xDir(1) == -1 && xDir(2) == 0 && xDir(3) == 0
    matRad_cfg.dispInfo('\nMirroring x-direction...')
    origCt = flip(origCt,1);
    matRad_cfg.dispInfo('finished!\n')
else
    nonStandardDirection = true;
end
    
if yDir(1) == 0 && yDir(2) == 1 && yDir(3) == 0
    matRad_cfg.dispInfo('y-direction OK\n')
elseif yDir(1) == 0 && yDir(2) == -1 && yDir(3) == 0
    matRad_cfg.dispInfo('\nMirroring y-direction...')
    origCt = flip(origCt,2);
    matRad_cfg.dispInfo('finished!\n')
else
    nonStandardDirection = true;
end

if nonStandardDirection
    matRad_cfg.dispWarning(['Non-standard patient orientation.\n'...
        'CT might not fit to contoured structures\n'])
end

%% interpolate cube
matRad_cfg.dispInfo('\nInterpolating CT cube...');
if ~isempty(obj.ImportGrid)
    obj.ct = matRad_interpDicomCtCube(origCt, obj.importCT.ctInfo, obj.importCT.resolution, obj.ImportGrid);
else
    obj.ct = matRad_interpDicomCtCube(origCt, obj.importCT.ctInfo, obj.importCT.resolution);
end
matRad_cfg.dispInfo('finished!\n');

%% remember some parameters of original dicom
tmp = [obj.importCT.ctInfo.ImagePositionPatient];

obj.ct.dicomInfo.PixelSpacing            = obj.importCT.ctInfo(1).PixelSpacing;                                      
obj.ct.dicomInfo.SlicePositions          = tmp(3,:);
obj.ct.dicomInfo.SliceThickness          = obj.importCT.ctInfo(1).SliceThickness;
obj.ct.dicomInfo.ImagePositionPatient    = obj.importCT.ctInfo(1).ImagePositionPatient;
obj.ct.dicomInfo.ImageOrientationPatient = obj.importCT.ctInfo(1).ImageOrientationPatient;
obj.ct.dicomInfo.PatientPosition         = obj.importCT.ctInfo(1).PatientPosition;
obj.ct.dicomInfo.Width                   = obj.importCT.ctInfo(1).Width;
obj.ct.dicomInfo.Height                  = obj.importCT.ctInfo(1).Height;
obj.ct.dicomInfo.RescaleSlope            = obj.importCT.ctInfo(1).RescaleSlope;
obj.ct.dicomInfo.RescaleIntercept        = obj.importCT.ctInfo(1).RescaleIntercept;
if isfield(completeDicom, 'Manufacturer')
obj.ct.dicomInfo.Manufacturer            = completeDicom.Manufacturer;
end
if isfield(completeDicom, 'ManufacturerModelName')
obj.ct.dicomInfo.ManufacturerModelName   = completeDicom.ManufacturerModelName;
end
if isfield(completeDicom, 'ConvolutionKernel')
obj.ct.dicomInfo.ConvolutionKernel       = completeDicom.ConvolutionKernel;
end

% store patientName only if user wants to
if isfield(completeDicom,'PatientName') && obj.dicomMetaBool == true
    obj.ct.dicomInfo.PatientName         = completeDicom.PatientName;
end
if obj.dicomMetaBool == true
    obj.ct.dicomMeta                     = completeDicom;
end

obj.ct.timeStamp = datestr(clock);

% convert to Hounsfield units
matRad_cfg.dispInfo('\nconversion of ct-Cube to Hounsfield units...');
obj = matRad_calcHU(obj);
matRad_cfg.dispInfo('finished!\n');

end
