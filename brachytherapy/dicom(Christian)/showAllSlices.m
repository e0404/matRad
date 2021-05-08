function showAllSlices(obj,modality,item)
%SHOWALLSLICES Summary of this function goes here
%   Detailed explanation goes here

noSlices = length(obj.(modality).(item).z);



for i=1:noSlices
    showSlice(obj,item,modality,i,'px','transversal');
    title(['transversal slice (',num2str(obj.(modality).(item).z(i),'%4.2f'),' [mm])'])
    pause(0.01)
    drawnow
end


end

