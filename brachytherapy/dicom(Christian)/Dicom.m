classdef Dicom < handle
    %DICOM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        CT
        US
        MR
        RS
        RP
        additionalInfo
        
        contourInterpolationRes
        
    end
    
    properties (SetAccess=private)
        CTfiles
        USfiles
        MRfiles
        RSfiles
        RPfiles
        changedFlag
        
        erasedRS;
    end
    
    methods
        function obj = Dicom(foldername,useOnly)
            
            global mm;
            
            if nargin<2
                useOnly = [];
            end
            
            if isdir(foldername)==0
                error('Folder does not exist')
            end
            
            obj.contourInterpolationRes = 1*mm;
           
            
            


            files = dir(foldername);
            files = files(~(strcmp({files(:).name},'.')' | strcmp({files(:).name},'..')'));
            
            testDir = struct2cell(files);
            testDir = cell2mat(testDir(4,1:end));
            
            
            if all([files(:).isdir])
                for i=1:length(files)
                    actFolder = [foldername,'/',files(i).name];
                    if i==1;
                        obj.analyzeFolder(actFolder,useOnly);
                    else
                        obj.addNewFolder(actFolder);
                    end
                end
            else
                obj.analyzeFolder(foldername,useOnly);
            end
            
            if ~isempty(obj.CTfiles) && isempty(useOnly)
                obj.buildVolumes('CT');
            elseif ~isempty(obj.CTfiles) && ~isempty(useOnly)
                if any(cellfun(@(x) strcmp(x,'CT'),useOnly))
                    obj.buildVolumes('CT');
                else
                    obj.changedFlag(1) = false;
                end
            end
            
            if ~isempty(obj.USfiles) && isempty(useOnly)
                obj.buildVolumes('US');
            elseif ~isempty(obj.USfiles) && ~isempty(useOnly)
                if any(cellfun(@(x) strcmp(x,'US'),useOnly))
                    obj.buildVolumes('US');
                end
            end
            
            if ~isempty(obj.MRfiles) && isempty(useOnly)
                obj.buildVolumes('MR');
            elseif ~isempty(obj.MRfiles) && ~isempty(useOnly)
                if any(cellfun(@(x) strcmp(x,'MR'),useOnly))
                    obj.buildVolumes('MR');
                end
            end
            
            if ~isempty(obj.RPfiles) && isempty(useOnly)
                obj.extractPlans();
            elseif ~isempty(obj.RPfiles) && ~isempty(useOnly)
                if any(cellfun(@(x) strcmp(x,'RP'),useOnly))
                    obj.extractPlans();
                end
            end
            
            if ~isempty(obj.RSfiles) && isempty(useOnly)
                obj.buildVOIs();
            elseif ~isempty(obj.RSfiles) && ~isempty(useOnly)
                if any(cellfun(@(x) strcmp(x,'RS'),useOnly))
                    obj.buildVOIs();
                end
            end
        end
        
        [mx,my,mz] = showCTSlice(obj,item,px,unit,orientation)
        [mx,my,mz] = showUSSlice(obj,item,px,unit,orientation)
        [mx,my,mz] = showSlice(obj,item,modality,px,unit,orientation)
        [imx,imy,imz,dose] = showIsodoseInSlice(obj,item,px,orientation,isodose,plan,source)
        [legendStr,contour] = showContourInSlice(obj,item,pos,orientation,showPlot)
        
        
        
        showAllSlices(obj,modality,item);
        
        txt = show3DStructureSets(obj,item,organ)
        
        function discardROI(obj,rsItem,roi)
            if isfield(obj.RS.(rsItem).structures,roi)
                obj.erasedRS.(roi) = obj.RS.(rsItem).structures.(roi);
                obj.RS.(rsItem).structures = rmfield(obj.RS.(rsItem).structures,roi);
                fprintf(' * %s discarded\n',roi)
            end
        end
        
        function addDose(obj,x,y,z,dose)
            if isempty(obj.RD)
                actItem = 'Item_1';
            else
                actItem = ['Item_',length(fieldnames(obj.RD))];
                obj.RD.(actItem).x = x;
                obj.RD.(actItem).y = y;
                obj.RD.(actItem).z = z;
                obj.RD.(actItem).dose = dose;
            end
        end
        
        function addNewFolder(obj,foldername)
            obj.analyzeFolder(foldername);
            
            if  obj.changedFlag(1) == true;
                obj.CT = [];
                obj.buildVolumes('CT');
            end
            
            if  obj.changedFlag(2) == true;
                obj.RS = [];
                obj.buildVOIs();
            end
            
            if  obj.changedFlag(3) == true;
                obj.US = [];
                obj.buildVolumes('US');
            end
            
            if  obj.changedFlag(4) == true;
                obj.MR = [];
                obj.buildVolumes('MR');
            end
            
            if  obj.changedFlag(5) == true;
                obj.RP = [];
                obj.extractPlans();
            end
        end
        
        
        function needleID = getNeedleMat(obj,item)
            needleID = zeros(length(obj.RP.(item).seeds.dwellPosition),3);
            actInd = 1;
            for i=1:length(obj.RP.(item).seeds.dwellPosition)
                
                needleID(i,:) = [actInd,actInd+size(obj.RP.(item).seeds.dwellPosition{i},1)-1,i];
                
                actInd = actInd+size(obj.RP.(item).seeds.dwellPosition{i},1);
            end
        end
        
        function [volItem,rsItem] = getVolumeAndStructures(obj,item)
            
            volItem=[];
            rsItem = [];
            imgID  = [];
            rsNames =  fieldnames(obj.RS);
            irs=0;
            found = 0;
            while irs<length(rsNames) && ~found
                irs = irs+1;
                found = strcmp(obj.RS.(rsNames{irs}).additionalInfo.structID,obj.RP.(item).additionalInfo.structID);
            end
            if found
                rsItem = rsNames{irs};
                imgID  = obj.RS.(rsNames{irs}).additionalInfo.imgID;
            end
            
            if ~isempty(imgID)
                modalities = {'US','CT','MR'};
                im         = 0;
                found      = 0;
                
                while im<3 && ~found
                    im = im+1;
                    if ~isempty(obj.(modalities{im}))
                        vNames = fieldnames(obj.(modalities{im}));
                        iv     = 0;
                        while iv<length(vNames) && ~found
                            iv = iv+1;
                            found = strcmp(obj.(modalities{im}).(vNames{iv}).additionalInfo.imgID,imgID);
                        end
                    end
                end
                if found
                    volItem{1} = modalities{im};
                    volItem{2} = vNames{iv};
                end
            end
        end
        
        function contour = getPrecalculatedContours(obj,rpItem)
            
            [volItem,rsItem] = obj.getVolumeAndStructures(rpItem);
            
            
            coordinates.x = obj.(volItem{1}).(volItem{2}).x;
            coordinates.y = obj.(volItem{1}).(volItem{2}).y;
            coordinates.z = obj.(volItem{1}).(volItem{2}).z;
            
            sNames = fieldnames(obj.RS.(rsItem).structures);
            contour = struct('maps',[],'color',[],'name',[]);
            for is = 1:length(sNames)
                if isempty(obj.RS.(rsItem).structures.(sNames{is}).additionalContours)
                    obj.RS.(rsItem).structures.(sNames{is}).precalcContours(coordinates);
                end
                contour(is).maps          = obj.RS.(rsItem).structures.(sNames{is}).additionalContours;
                contour(is).color         = obj.RS.(rsItem).structures.(sNames{is}).color;
                contour(is).name          = sNames{is};
            end
        end
        
        function [x,y,z] = getObjectCoordinates(obj,name,rpItem,margin)
            
            if nargin<4
                margin = 0;
            end
            
            [volItem,rsItem] = obj.getVolumeAndStructures(rpItem);
            
            
            x = obj.(volItem{1}).(volItem{2}).x;
            y = obj.(volItem{1}).(volItem{2}).y;
            z = obj.(volItem{1}).(volItem{2}).z;
            
            sNames = fieldnames(obj.RS.(rsItem).structures);
            found = false;
            is    = 0;
            while is<length(sNames) && ~found
                is    = is+1;
                found = strcmp(name,obj.RS.(rsItem).structures.(sNames{is}).name);
            end
            x = x(obj.RS.(rsItem).structures.(sNames{is}).minXYZ(1)-margin<=x & x<=obj.RS.(rsItem).structures.(sNames{is}).maxXYZ(1)+margin);
            y = y(obj.RS.(rsItem).structures.(sNames{is}).minXYZ(2)-margin<=y & y<=obj.RS.(rsItem).structures.(sNames{is}).maxXYZ(2)+margin);
            z = z(obj.RS.(rsItem).structures.(sNames{is}).minXYZ(3)-margin<=z & z<=obj.RS.(rsItem).structures.(sNames{is}).maxXYZ(3)+margin);
            
        end
        
        function showCatheters(obj,rpItem,catNo)
            hold on
            if nargin<3
                for in = 1:length(obj.RP.(rpItem).seeds.dwellPosition)
                    plot3(obj.RP.(rpItem).seeds.dwellPosition{in}(:,1),obj.RP.(rpItem).seeds.dwellPosition{in}(:,2),obj.RP.(rpItem).seeds.dwellPosition{in}(:,3),'-o','Linewidth',2)
                end
            else
                plot3(obj.RP.(rpItem).seeds.dwellPosition{catNo}(:,1),obj.RP.(rpItem).seeds.dwellPosition{catNo}(:,2),obj.RP.(rpItem).seeds.dwellPosition{catNo}(:,3),'-o','Linewidth',2)
            end
        end
        
        extractPlans(obj)
        
        [dosePoints,volume,patch] = generateNormalTissue(obj,rpItem,additionalInfo)
        
        showContour(obj,rsItem,structName)
        
        posAtSlice = getNeedlePosition(obj,rpItem,posZ,noNeedle)
        
        [needles] = buildNeedles(obj,rpItem,needleDiameter)
        
        function fieldname = getStructuresFromVolume(obj,item,modality)
            sNames = fieldnames(obj.RS);
            is    = 0;
            found = false;
            while is<length(sNames) && ~found
                is = is+1;
                found = strcmp(obj.RS.(sNames{is}).additionalInfo.imgID,obj.(modality).(item).additionalInfo.imgID);
            end
            
            fieldname = sNames{is};
        end
        
        
        function item  = findPlan(obj,status)
            pNames = fieldnames(obj.RP);
            irp = 0;
            found = false;
            while irp <length(pNames) && ~found
                irp = irp+1;
                if ~isempty(obj.RP.(pNames{irp}).additionalInfo.status)
                    found = strcmpi(obj.RP.(pNames{irp}).additionalInfo.status,status);
                end
            end
            
            item = pNames{irp};
        end
        
        
    end
    
    methods (Access=private)
        analyzeFolder(obj,foldername,useOnly)
        buildCTVolumes(obj)
        buildUSVolumes(obj)
        buildVolumes(obj,modality)
        buildVOIs(obj)
        
    end
    
    methods (Static)
        [contPoints,contInd,splittedContPoints] = intersectPlanObject(planeNormal, planeCenter, objectVertices, objectFaces)
        [cluster,z,logVolume] = clusterData(data,dz)
    end
    
end

