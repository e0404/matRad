function [needles] = buildNeedles(obj,rpItem,needleDiameter)
%BUILDNEEDLES Summary of this function goes here
%   Detailed explanation goes here

global mm

noPoints = 10;
theta = linspace(0,2*pi,noPoints+1)';
rho   = ones(size(theta))*needleDiameter/2;
[x,y] = pol2cart(theta(1:end-1),rho(1:end-1));

if isfield(obj.RP.(rpItem),'needles')
    disp('* Needles found')
    noNeedles =length(obj.RP.(rpItem).needles.pos);
    needles = cell(noNeedles,1);

        
    %%
    for in = 1:noNeedles
        pos   = obj.RP.(rpItem).needles.pos{in};
        tempZ = pos(:,3);
        
                 map = cell(length(tempZ),1);

            for ip=1:length(map)
                 map{ip} = [pos(ip,1)+x,pos(ip,2)+y,ones(size(x))*pos(ip,3)];
            end

            
            outerContourMap.contourMap = map;
            outerContourMap.z = tempZ;
            volume          = [];
            innerContourMap = [];
            
            %%
            fv.faces = [];
            fv.vertices = [];
            shift = 0;
            for i=1:length(map)-1
                pos = [map{i};map{i+1}];
                fv.faces  = [fv.faces;shift+[1:noPoints-1;noPoints+1:2*noPoints-1;noPoints+2:2*noPoints]';shift+[1:noPoints-1;noPoints+2:2*noPoints;2:noPoints]'];
                
                
                
                fv.vertices  = [fv.vertices;pos];
                shift = size(fv.vertices,1);
            end
            
            
            %[outerContourMap,innerContourMap,volume,fv] = analyzeContours(tempZ,map);
            fv.FaceColor = [0.4 0.4 0.4];
            fv.FaceAlpha = 0.3;
            fv.EdgeColor = 'none';
            needles{in}  = newVOI('needle',fv,outerContourMap,innerContourMap,volume);
    end
else
    
    
    volItem = getVolumeAndStructures(obj,rpItem);
    noNeedles = length(obj.RP.(rpItem).seeds.dwellPosition);
    
    if isempty(volItem{1})
        %%
        sNames = fieldnames(obj.RS.Item_1.structures);
        zMinMax = [inf,-inf];
        for i=1:length(sNames)
            if obj.RS.Item_1.structures.(sNames{i}).minXYZ(3)<zMinMax(1)
                zMinMax(1) = obj.RS.Item_1.structures.(sNames{i}).minXYZ(3);
            end
            if obj.RS.Item_1.structures.(sNames{i}).maxXYZ(3)>zMinMax(2)
                zMinMax(2) = obj.RS.Item_1.structures.(sNames{i}).maxXYZ(3);
            end
        end
        
        z = linspace(zMinMax(1),zMinMax(2),32)';
        
    else
        z = obj.(volItem{1}).(volItem{2}).z;
    end
    
    
    needles = cell(noNeedles,1);
    foundNeedles = 0;
    for in=1:noNeedles;
        if ~isempty(obj.RP.(rpItem).seeds.dwellPosition{in}) && size(obj.RP.(rpItem).seeds.dwellPosition{in},1)>1
            minZ = min(obj.RP.(rpItem).seeds.dwellPosition{in}(:,3))-obj.RP.(rpItem).additionalInfo.stepSize(in)/2;
            maxZ = max(obj.RP.(rpItem).seeds.dwellPosition{in}(:,3))+obj.RP.(rpItem).additionalInfo.stepSize(in)/2;
            
            tempZ = z(minZ<=z & z<=maxZ);
            pos = getPointsAtNeedle(obj.RP.(rpItem).seeds.dwellPosition{in},[],tempZ);
            %%pos = getNeedlePosition(obj,rpItem,tempZ,in);
            
            map = cell(length(pos),1);
            for ip=1:length(pos)
                map{ip} = [pos(ip,1)+x,pos(ip,2)+y,ones(size(x))*pos(ip,3)];
            end
            
            
            [outerContourMap,innerContourMap,volume,fv] = analyzeContours(tempZ,map);
            fv.FaceColor = [0.4 0.4 0.4];
            fv.FaceAlpha = 0.3;
            fv.EdgeColor = 'none';
            foundNeedles = foundNeedles+1;
            needles{foundNeedles}  = newVOI('needle',fv,outerContourMap,innerContourMap,volume);
        else
            if isempty(obj.RP.(rpItem).seeds.dwellPosition{in})
                disp('Needle is empty')
            end
            
            %if size(obj.RP.(rpItem).seeds.dwellPosition{in},1)==1
            %    disp('Needle is ignored')
            %end
            
        end
        
    end
    
    needles = needles(1:foundNeedles);
end


end

