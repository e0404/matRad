function txt = show3DStructureSets(obj,item,organ)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
hold on
txt = [];
if nargin<3
    
    fnames = fieldnames(obj.RS.(item).structures);
    
    legendTxt = [];
    for i=1:length(fnames)
        legendTxt{i} = fnames{i};
        if strcmp(fnames{i},'ptv_high')
            legendTxt{i} = 'PTV';
            
            obj.RS.(item).structures.(fnames{i}).showObject('b');
        elseif strcmp(fnames{i},'ptv_low')
            legendTxt{i} = 'CTV';
            obj.RS.(item).structures.(fnames{i}).showObject('r');
        else
            obj.RS.(item).structures.(fnames{i}).showObject();
        end
        
        
    end
    obj.RS.(item).structures.(fnames{i}).standardPlotSetting;
    legend(legendTxt, 'Interpreter', 'none');
    txt = legendTxt;
else
    
    obj.RS.(item).structures.(organ).showObject();
    obj.RS.(item).structures.(organ).standardPlotSetting;
end
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
end

