function [imx,imy,imz,dose] = showIsodoseInSlice(obj,item,px,orientation,isodose,plan,source)
%SHOWISODOSEINSLICE Summary of this function goes here
%   Detailed explanation goes here

hold on
if ~isstruct(plan)
    dose = plan;
else
    dose = [];
end
if isempty(obj.CT)  
    modality = 'US';
elseif isempty(obj.US)
    modality = 'CT';
else
    error('Modality needs to be specified')
end

c = jet(length(isodose));
 
switch orientation
    case 'transversal'
        
        if isempty(dose)
        [imx,imy] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y);
         imz = ones(size(imx(:,:,1)))*px;
        
        
%         slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,mx(:,:,1),my(:,:,1),ones(size(mx(:,:,1)))*px);
     
dose(:,:)  = calculateDoseInVolume(imx,imy,imz,plan,source);

        
        end
        
%         imagesc(obj.(modality).(item).x,obj.(modality).(item).y,slice,[-200 200]);
%         colormap('gray')
%         hold on

       
        for id=1:length(isodose)
            contour(obj.(modality).(item).x,obj.(modality).(item).y,dose,[isodose(id) isodose(id)],'LineColor',c(id,:),'LineWidth',1.5)
        end
        
    case 'frontal'
        
        
        %[mx,my,mz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y,obj.(modality).(item).z);
        if isempty(dose)
        [imx,imz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).z);
        imy = ones(size(imz))*px;
        %slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,imx,ones(size(imz))*px,imz);
        dose(:,:)  = calculateDoseInVolume(imx,imy,imz,plan,source);
        
        end
        
        
%         imagesc(obj.(modality).(item).x,obj.(modality).(item).z,slice,[-200 200]);
%         colormap gray
%         hold on
%         
    
        for id=1:length(isodose)
            contour(obj.(modality).(item).x,obj.(modality).(item).z,dose,[isodose(id) isodose(id)],'LineColor',c(id,:))
        end
        
%         hold off
%         set(gca,'YDir','normal');
%         
%         axis image
        
    case 'sagittal'
        
        %[mx,my,mz] = meshgrid(obj.(modality).(item).x,obj.(modality).(item).y,obj.(modality).(item).z);
         if isempty(dose)
        [imy,imz] = meshgrid(obj.(modality).(item).y,obj.(modality).(item).z);
        
        imx = ones(size(imy))*px;
        %slice(:,:) = interp3(mx,my,mz,obj.(modality).(item).volume,ones(size(imz))*px,imy,imz);
        dose(:,:)  = calculateDoseInVolume(imx,imy,imz,plan,source);

         end
        
%         imagesc(obj.(modality).(item).y,obj.(modality).(item).z,slice,[-200 200]);
%         hold on

        for id=1:length(isodose)
            contour(obj.(modality).(item).y,obj.(modality).(item).z,dose,[isodose(id) isodose(id)],'LineColor',c(id,:))
        end
%         hold off
%         colormap gray
%         axis image
%         set(gca,'YDir','normal');
    otherwise
        error('Orientation not specified: transversal, frontal, sagittal');
        
end

hold off
end


% function [dose] = calculateDoseInVolume(mx,my,mz,plan,source)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% 
% dose = zeros(size(mx));
% 
% for ip=1:length(plan.dwellTime)
%     source.setPosition(plan.dwellPosition(ip,:));
%     source.setDirection(plan.dwellDirection(ip,:));
%     
%     dose(:) = dose(:)+plan.dwellTime(ip)*source.getDoseDictionary(mx(:),my(:),mz(:));
% end
% 
% 
% 
% end