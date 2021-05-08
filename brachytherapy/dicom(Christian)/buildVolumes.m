function buildVolumes(obj,modality)
%BUILDUSVOLUMES Summary of this function goes here
%   Detailed explanation goes here

%%
%prec  = 1e-9;
files = [modality,'files'];

% if ~isempty(obj.(modality))
%     itemOffset = length(fieldnames(obj.(modality)));
% else
%     itemOffset = 0;
% end

%seriesNumber            = zeros(length(obj.(files)),1);
%instanceNumber          = zeros(length(obj.(files)),1);
imagePositionPatient    = zeros(length(obj.(files)),3);
imageOrientationPatient = zeros(length(obj.(files)),6);
pixelSpacing            = zeros(length(obj.(files)),2);
rows                    = zeros(length(obj.(files)),1);
columns                 = zeros(length(obj.(files)),1);
seriesInst              = cell(length(obj.(files)),1);
volStruct               = cell(length(obj.(files)),1);
frameOfReference        = cell(length(obj.(files)),1);
h = waitbar(0,'Please wait processing image data...');
steps = length(obj.(files));
actImg = 0;
for i=1:length(obj.(files))
    inf                           = dicominfo(obj.(files){i});
    if isfield(inf,'ImagePositionPatient')
      
        actImg  = actImg +1;
        
    %seriesNumber(actImg)              = inf.SeriesNumber;
    %instanceNumber(actImg)            = inf.InstanceNumber;
    imagePositionPatient(actImg,:)     = inf.ImagePositionPatient;
    imageOrientationPatient(actImg,:)  = inf.ImageOrientationPatient;
    seriesInst{actImg}                 = inf.SeriesInstanceUID;
    frameOfReference{actImg}           = inf.FrameOfReferenceUID;
    pixelSpacing(actImg,:)             = inf.PixelSpacing;
    rows(actImg)                      = inf.Rows;
    columns(actImg)                   = inf.Columns;
    
    
    slice = dicomread(obj.(files){actImg});
    if strcmp(modality,'CT')
        slice = double(inf.RescaleSlope)*double(slice)+double(inf.RescaleIntercept);
    end
    volStruct{actImg} = slice;

    end
        waitbar(i / steps)
end
close(h)

%%
%seriesNumber            = seriesNumber(1:actImg) ;
%instanceNumber          = instanceNumber(1:actImg) ;
imagePositionPatient    = imagePositionPatient(1:actImg,:);
imageOrientationPatient = imageOrientationPatient(1:actImg,:);
seriesInst              = seriesInst(1:actImg);
frameOfReference        = frameOfReference(1:actImg);
pixelSpacing            = pixelSpacing(1:actImg,:); 
rows                    = rows(1:actImg);
columns                 = columns(1:actImg);
volStruct               =  volStruct(1:actImg);

%%
sID = unique(seriesInst);

for isd=1:length(sID)
    % get all images with the same ID
    ind = find(ismember(seriesInst,sID{isd}));
    %extract position
    actPos = imagePositionPatient(ind,:);
    
    % check if all images have the same dimesnion and continue
    if ~(all(columns(ind)==columns(ind(1))) && all(rows(ind)==rows(ind(1))))
        error('Images should have the same dimension')
    end 
    actRows   = rows(ind(1));
    actColums = columns(ind(1));
    actSlices = length(ind);
    
    % save images to volume
    volume = zeros(actRows,actColums,actSlices);
    for i=1:actSlices 
         volume (:,:,i) = volStruct{ind(i)};
    end
    
    % extract the coordinates
    tempX = zeros(actSlices,actColums);
    tempY = zeros(actSlices,actRows);
    for i=1:actSlices  
        M = zeros(3,4);
        M(:,1:3) = [imageOrientationPatient(ind(i),1:3)'*pixelSpacing(ind(i),1),imageOrientationPatient(ind(i),4:6)'*pixelSpacing(ind(i),2),zeros(3,1)];
        M(:,4) = imagePositionPatient(ind(i),:)';
        ix = 0:actColums-1;
        jy = 0:actRows-1;
        
        % if the coordinate system is not orthorgonal throw an exception
        if ~(all(M(1,2:3)==0) &&  all(M(2,[1,3])==0) && all(M(3,1:3)==0))
            disp('None orthogonal coordinate system')
        end
        
        % transform coordinates
        tempX(i,:) = M(1,1)*ix+M(1,4);
        tempY(i,:) = M(2,2)*jy+M(2,4);
    end
    
    % get coordinate system and check if it is valid
    ux        = unique(tempX,'rows');
    uy        = unique(tempY,'rows');
    [uz,indZ] = sort(actPos(:,3));
    
    if size(ux,1)>1 || size(uy,1)>1
        error('Coordinates vary for different slices')
    end
    
    % saving the volume
    obj.(modality).(['Item_',num2str(isd)]).volume = volume(:,:,indZ);
    obj.(modality).(['Item_',num2str(isd)]).x = ux(:);
    obj.(modality).(['Item_',num2str(isd)]).y = uy(:);
    obj.(modality).(['Item_',num2str(isd)]).z = uz(:);
    obj.(modality).(['Item_',num2str(isd)]).additionalInfo.imgID    = frameOfReference{ind(1)};
    obj.(modality).(['Item_',num2str(isd)]).additionalInfo.seriesID = sID{isd};
    obj.(modality).(['Item_',num2str(isd)]).additionalInfo.filename = obj.(files)(ind);
end
    
   
% if length(unique(imagePositionPatient(:,1)))==1 && length(unique(imagePositionPatient(:,2)))==1
%     
%     
%     z                    = imagePositionPatient(:,3);
%     imagePositionPatient = unique(imagePositionPatient(:,1:2),'rows');
%     columns              = unique(columns);
%     rows                 = unique(rows);
%     
%     pixelSpacing = unique(pixelSpacing,'rows');
%     imageOrientationPatient = unique(imageOrientationPatient,'rows');
%     
%     M = [imageOrientationPatient(1:3)'*pixelSpacing(1),imageOrientationPatient(4:6)'*pixelSpacing(2),zeros(3,1)];
%     M = [M(1:2,:),imagePositionPatient'];
%     i = 0:columns-1;
%     j = 0:rows-1;
%     if length(i)==length(j)
%         
%         iVec = [i;j;zeros(1,rows);ones(1,rows)];
%         vec = M*iVec;
%         if find(M(:,1))==1 && find(M(:,2))==2
%             x = vec(1,:)';
%             y = vec(2,:)';
%         else
%             error('Unknown axis configuration detected')
%         end
%         
%         
%     else
%         if abs(M(1,2))<prec && abs(M(2,1))<prec
%             x = (M(1,1)*i+M(1,4))';
%             y = (M(2,2)*j+M(2,4))';
%         else
%             disp('test')
%             error('Size of rows and columns do not match')
%         end
%         
%     end
%     
%     volume = zeros(rows,columns,length(obj.(files)));
%     for i=1:length(obj.(files))
%         volume (:,:,i) = volStruct{i};
%     end
%     uSI = unique(seriesInst);
%     if length(uSI)>1 || length(unique(seriesNumber))>1
%         for i=1:length(uSI)
%             indID = ismember(seriesInst,uSI{i});
%             [tempz,ind] = sort(z(indID));
%             obj.(modality).(['Item_',num2str(i)]).volume = volume(:,:,indID);
%             obj.(modality).(['Item_',num2str(i)]).volume = obj.(modality).(['Item_',num2str(i)]).volume(:,:,ind);
%             obj.(modality).(['Item_',num2str(i)]).x = x;
%             obj.(modality).(['Item_',num2str(i)]).y = y;
%             obj.(modality).(['Item_',num2str(i)]).z                 = tempz;
%             obj.(modality).(['Item_',num2str(i)]).additionalInfo.imgID = inf.FrameOfReferenceUID;
%         end
%     else
%         [z,ind] = sort(z);
%         obj.(modality).Item_1.volume = volume(:,:,ind);
%         obj.(modality).Item_1.x      = x;
%         obj.(modality).Item_1.y      = y;
%         obj.(modality).Item_1.z      = z;
%         obj.(modality).Item_1.additionalInfo.imgID = inf.FrameOfReferenceUID;
%     end
% else
%     
%     uniqueSeriesInst = unique(seriesInst);
%     %%
%     for iu = 1:size(uniqueSeriesInst,1)
%         disp('Image position patient is not unique')
%         indSeries = find(ismember(seriesInst,uniqueSeriesInst{1}));
%         z  = imagePositionPatient(indSeries,3);
%         [z,ind] = sort(z);
%         
%         volume = zeros(rows(indSeries(1)),columns(indSeries(1)),length(z));
%         for i=1:length(z)
%             volume (:,:,i) = volStruct{indSeries(i)};
%         end
%         
%         tempPixelSpacing = unique(pixelSpacing(indSeries,:),'rows');
%         tempImageOrientationPatient = unique(imageOrientationPatient(indSeries,:),'rows');
%         
%         M = [tempImageOrientationPatient(1:3)'*tempPixelSpacing(1),tempImageOrientationPatient(4:6)'*tempPixelSpacing(2),zeros(3,1)];
%         for iz=1:length(z)
%             M = [M(1:2,1:3),imagePositionPatient(indSeries(iz),1:2)'];
%             i = 0:columns(indSeries(iz))-1;
%             j = 0:rows(indSeries(iz))-1;
%             if length(i)==length(j)
%                 iVec = [i;j;zeros(1,rows(indSeries(iz)));ones(1,rows(indSeries(iz)))];
%                 vec = M*iVec;
%                 x   = vec(1,:)';
%                 y   = vec(2,:)';
%             else
%                 if abs(M(1,2))<prec && abs(M(2,1))<prec
%                     x = (M(1,1)*i+M(1,4))';
%                     y = (M(2,2)*j+M(2,4))';
%                 else
%                     disp('test')
%                     error('Size of rows and columns do not match')
%                 end
%                 
%             end
%             obj.(modality).(['Item_',num2str(iu+itemOffset)]).x{iz} = x;
%             obj.(modality).(['Item_',num2str(iu+itemOffset)]).y{iz} = y;
%             
%         end
%         
%         obj.(modality).(['Item_',num2str(iu+itemOffset)]).z      = z;
%         obj.(modality).(['Item_',num2str(iu+itemOffset)]).volume =  volume(:,:,ind);
%         obj.(modality).(['Item_',num2str(iu+itemOffset)]).additionalInfo.imgID = frameOfReference{indSeries(1)};
%         
%         tempX = unique(cell2mat(obj.(modality).(['Item_',num2str(iu+itemOffset)]).x)','rows');
%         tempY = unique(cell2mat(obj.(modality).(['Item_',num2str(iu+itemOffset)]).y)','rows');
%         if size(tempX,1) ==1
%             obj.(modality).(['Item_',num2str(iu+itemOffset)]).x = tempX(:);
%         else
%             disp('x not unique')
%         end
%         if size(tempY,1) ==1
%             obj.(modality).(['Item_',num2str(iu+itemOffset)]).y = tempY(:);
%         else
%             disp('y not unique')
%         end
%     end
% end

end

