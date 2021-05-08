function obj = buildVOIs( obj )
%BUILDVOIS Summary of this function goes here
%   Detailed explanation goes here


for is=1:length(obj.RSfiles)
    
    fprintf('Object %i of %i\n',is,length(obj.RSfiles))
    
    inf = dicominfo(obj.RSfiles{is});
    itemName = ['Item_',num2str(is)];
    
    sNames = fieldnames(inf.StructureSetROISequence);
    
    h = waitbar(0,'Please wait processing VOI data...');
    steps = length(sNames);
    for sn=1:length(sNames)
        actI = sNames{sn};
        name  = inf.StructureSetROISequence.(actI).ROIName;
        name = lower(strrep(name,' ','_'));
        name = strrep(name,'.','');
        name = strrep(name,'*','');
        
        if isfield(inf.ROIContourSequence.(actI),'ContourSequence')
            
            if ~isfield(inf.ROIContourSequence.(actI),'ROIDisplayColor')
                color = [1 0 0];
            else
                color = inf.ROIContourSequence.(actI).ROIDisplayColor;
            end
            
            % get contour map
            fNames = fieldnames(inf.ROIContourSequence.(actI).ContourSequence);
            
            map = cell(length(fNames),1);
            z = zeros(length(fNames),1);
            
            for i=1:length(fNames)
                
                icn =i;% inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).ContourNumber changed on 05.08.2015
                icp = inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).NumberOfContourPoints;
                size(inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).ContourData);
                idX = 1:3:3*icp;
                idY = idX+1;
                idZ = idX+2;
                
                
                
                map{icn} = [inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).ContourData(idX),...
                    inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).ContourData(idY),...
                    inf.ROIContourSequence.(actI).ContourSequence.(fNames{i}).ContourData(idZ)];
                
                if size(map{icn},1)~=icp
                    error('Points are missing')
                end
                
                z(icn) = map{icn}(1,3);
            end
            
            if length(z)>1
                if ~isempty(obj.contourInterpolationRes)
                    [outerContourMap,innerContourMap,volume,fv] = analyzeContours(z,map,obj.contourInterpolationRes);
                else
                    [outerContourMap,innerContourMap,volume,fv] = analyzeContours(z,map);
                end
                
                
                fv.FaceColor = color./256;
                fv.FaceAlpha = 0.3;
                fv.EdgeColor = 'none';
                
                name = strrep(name,'-','_');
                name = strrep(name,'+','_');
                name = strrep(name,'#','_');
                obj.RS.(itemName).structures.(name) = newVOI(name,fv,outerContourMap,innerContourMap,volume);
                obj.RS.(itemName).additionalInfo.structID = inf.SOPInstanceUID;
                obj.RS.(itemName).additionalInfo.imgID    = inf.ReferencedFrameOfReferenceSequence.Item_1.FrameOfReferenceUID;
                obj.RS.(itemName).additionalInfo.filename = obj.RSfiles{is};
            else
                fprintf('ROI with only one contour\n');
            end
        else
            fprintf('Unable to find contour sequence in ROI: %s\n',name);
        end
        waitbar(sn / steps)
    end
    delete(h)
    
end
%obj = Object3D(name,faces,vertices,color,alpha)

end



