function ct = matRad_parseMVF(ct,origResolution,fname,dir)

if nargin < 3
    [fname,dir] = uigetfile('.txt','Select image deformation vector file');
    if nargin < 2
        origResolution.x = str2double('What is the original resolution in the x direction?','s');
        origResolution.y = str2double('What is the original resolution in the y direction?','s');
        origResolution.z = str2double('What is the original resolution in the z direction?','s');
    end
end

oldDir = pwd;
cd(dir);

fileID = fopen(fname,'r');

%first line is empty
tline = fgetl(fileID);
%second line contains frame number
tline = fgetl(fileID);
spacesInd = strfind(tline,' ');
frame = str2double(tline(spacesInd(7):end));
%third and subsequent lines contain vector information
tline = fgetl(fileID);

%initialize motion vectors
ct.motionVecX{frame} = zeros(ct.cubeDim);
ct.motionVecY{frame} = zeros(ct.cubeDim);
ct.motionVecZ{frame} = zeros(ct.cubeDim);

%assume no motion for unspecified voxels
for i = 1:ct.cubeDim(1)
    ct.motionVecY{frame}(i,:,:) = i;
end
for i = 1:ct.cubeDim(2)
    ct.motionVecX{frame}(:,i,:) = i;
end
for i = 1:ct.cubeDim(3)
    ct.motionVecZ{frame}(:,:,i) = i;
end



resRatioX = ct.resolution.x./origResolution.x;
resRatioY = ct.resolution.y./origResolution.y;
resRatioZ = ct.resolution.z./origResolution.z;

if ~( round(resRatioX) == resRatioX && round(resRatioY) == resRatioY && round(resRatioZ) == resRatioZ )
    error('Non-integral ratios of resolutions.  Inexact motion vector mapping may occur.');
end

i = 1;

while ischar(tline)

    %{
    if ~strcmp(tline(1:5),'known')
        error('Unknown vector encountered');
    end
    %}
    
    spacesInd = strfind(tline,' ');
    
    % each line contains the voxel coordinate (in the original resolution)
    % of each transformed voxel, and also the voxel sub-coordinates in the
    % transformed frame
    xCoord_vox = 1+str2double(tline((spacesInd(2)+1):(spacesInd(3)-1)))./resRatioX;
    yCoord_vox = 1+str2double(tline((spacesInd(3)+1):(spacesInd(4)-1)))./resRatioY;
    zCoord_vox = 1+str2double(tline((spacesInd(4)+1):(spacesInd(5)-1)))./resRatioZ;
    
    if round(xCoord_vox) == xCoord_vox && round(yCoord_vox) == yCoord_vox && round(zCoord_vox) == zCoord_vox
        %only store motion vectors for voxels which map exactly to a voxel
        %in the interpolated resolution (note that the voxel in the new
        %frame will probably not map exactly to a voxel)
        xCoordNewFrame_vox = 1+str2double(tline((spacesInd(7)+1):(spacesInd(8)-1)))./resRatioX;
        yCoordNewFrame_vox = 1+str2double(tline((spacesInd(8)+1):(spacesInd(9)-1)))./resRatioY;
        zCoordNewFrame_vox = 1+str2double(tline((spacesInd(9)+1):end))./resRatioZ;
        
        ct.motionVecX{frame}(yCoord_vox,xCoord_vox,zCoord_vox) = xCoordNewFrame_vox;
        ct.motionVecY{frame}(yCoord_vox,xCoord_vox,zCoord_vox) = yCoordNewFrame_vox;
        ct.motionVecZ{frame}(yCoord_vox,xCoord_vox,zCoord_vox) = zCoordNewFrame_vox;
        
        
        matRad_progress(i,numel(ct.motionVecX{frame}));
        i = i+1;
        tline = fgetl(fileID);
    else
        tline = fgetl(fileID);
        continue
    end
    
end

fclose(fileID);
cd(pwd);

end