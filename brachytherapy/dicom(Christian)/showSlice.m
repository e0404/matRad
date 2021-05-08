function [mx,my,mz] = showSlice(obj,item,modality,position,unit,orientation)
%SHOWSLICE Summary of this funUSion goes here
%   Detailed explanation goes here

prec = 1e-9;

switch nargin
    case 1
        item        = 'Item_1';
        position          = [];
        unit        = 'px';
        orientation = 'sagittal';
    case 2
        position          = [];
        unit        = 'px';
        orientation = 'transversal';
    case 3
        unit        = 'px';
        orientation = 'transversal';
    case 4
        orientation = 'transversal';
end


switch orientation
    case 'transversal'
        switch unit
            case 'px'
                if isempty(position)
                    position = round(size(obj.(modality).(item).volume,3)/2);
                end
                slice(:,:) = obj.(modality).(item).volume(:,:,position);
            case 'mm'
                
                
                dz = abs(obj.(modality).(item).z-position);
                ind = dz<prec;
                if any(ind)
                    slice(:,:) = obj.(modality).(item).volume(:,:,ind);
                else
                    [mx,my,mz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y,obj.(modality).(item).z);
                    slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,mx(:,:,1),my(:,:,1),ones(size(mx(:,:,1)))*position);
                end
        end
        
        imagesc(obj.(modality).(item).x,obj.(modality).(item).y,slice);
        xlabel('x [mm]')
        ylabel('y [mm]')
        colormap gray
        axis image
        
    case 'sagittal'
        
        switch unit
            case 'px'
                if isempty(position)
                    position = round(size(obj.(modality).(item).volume,3)/2);
                end
                slice(:,:) = obj.(modality).(item).volume(:,position,:);
                slice = fliplr(rot90(slice,-1));
            case 'mm'
                dy = abs(obj.(modality).(item).y-position);
                ind = dy<prec;
                if any(ind)
                    slice(:,:) = obj.(modality).(item).volume(:,position,:);
                    slice      = fliplr(rot90(slice,-1));
                else
                    
                    [mx,my,mz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y,obj.(modality).(item).z);
                    [imx,imz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).z);
                    
                    slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,imx,ones(size(imz))*position,imz);
                end
        end
        
        
        imagesc(obj.(modality).(item).x,obj.(modality).(item).z,slice);
        set(gca,'YDir','normal');
        xlabel('x [mm]')
        ylabel('z [mm]')
        colormap gray
        axis image
        
    case 'frontal'
        switch unit
            case 'px'
                if isempty(position)
                    position = round(size(obj.(modality).(item).volume,1)/2);
                end
                slice(:,:) = obj.(modality).(item).volume(position,:,:);
                slice = fliplr(rot90(slice,-1));
            case 'mm'
                
                dx = abs(obj.(modality).(item).x-position);
                ind = dx<prec;
                if any(ind)
                    slice(:,:) = obj.(modality).(item).volume(ind,:,:);
                    slice      = fliplr(rot90(slice,-1));
                else
                    [mx,my,mz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y,obj.(modality).(item).z);
                    
                    [imy,imz] = meshgrid(obj.(modality).(item).y,obj.(modality).(item).z);
                    
                    slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,ones(size(imz))*position,imy,imz);
                end
        end
        
        imagesc(obj.(modality).(item).y,obj.(modality).(item).z,slice);
        colormap gray
        axis image
        xlabel('y [mm]')
        ylabel('z [mm]')
        set(gca,'YDir','normal');
    otherwise
        error('Orientation not specified: transversal, frontal, sagittal');
        
end

