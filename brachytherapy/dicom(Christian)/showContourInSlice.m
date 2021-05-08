function [legendStr,contour] = showContourInSlice(obj,item,pos,orientation,showPlot)
%SHOWCONTOURINSLICE Summary of this function goes here
%   Detailed explanation goes here

if nargin<5
    showPlot = true;
end

switch orientation
    case 'transversal'
        position  = [0 0 pos];
        direction = [0 0 1];
    case 'frontal'
        position  = [0 pos 0];
        direction = [0 1 0];
    case 'sagittal'
        position  = [pos 0 0];
        direction = [1 0 0];
        
        
end
%clf
if showPlot
    hold on
end
iNames = fieldnames(obj.RS.(item).structures);
legendStr =  cell(length(iNames),1);
for in =1:length(iNames)
    [contPoints] = intersectPlaneObject(direction, position, obj.RS.(item).structures.(iNames{in}).patchStruct.vertices, obj.RS.(item).structures.(iNames{in}).patchStruct.faces);
    
    contour.(iNames{in}).contPoints = [];
    contour.(iNames{in}).color      = [];
    legendStr{in} = iNames{in};
    
    
    if ~isempty(contPoints)
        switch orientation
            case 'transversal'
                 ic=1;
                    if showPlot
                        h{ic} = plot(contPoints(:,1),contPoints(:,2),'Color',obj.RS.(item).structures.(iNames{in}).color);
                    end
                    contour.(iNames{in}).contPoints{ic} = [contPoints(:,1),contPoints(:,2)];
                    contour.(iNames{in}).color          = obj.RS.(item).structures.(iNames{in}).color;
                
                
            case 'frontal'
                ic=1;
                    if   showPlot
                        h{ic} = plot(contPoints(:,1),contPoints(:,3),'Color',obj.RS.(item).structures.(iNames{in}).color);
                    end
                    contour.(iNames{in}).contPoints{ic} = [contPoints(:,1),contPoints(:,3)];
                    contour.(iNames{in}).color          = obj.RS.(item).structures.(iNames{in}).color;
  
            case 'sagittal'
                ic=1;
                    if showPlot
                        h{ic} = plot(contPoints(:,2),contPoints(:,3),'Color',obj.RS.(item).structures.(iNames{in}).color);
                    end
                    contour.(iNames{in}).contPoints{ic} = [contPoints(:,2),contPoints(:,3)];
                    contour.(iNames{in}).color          = obj.RS.(item).structures.(iNames{in}).color;

        end
    end
    if showPlot

            set(h{1},'linewidth',3);
            set(h{1},'linestyle','--');

    end
end

hold off

end


