matRad_rc

load BOXPHANTOM

ct.cubeDim = [320,320,320];

lowX    =   1;
lowY    =   1;
lowZ    =   1;
highX   = 320;
highY   = 320;
highZ   = 320;

indices = zeros((highX-lowX+1)*(highY-lowY+1)*(highZ-lowZ+1),3);

position = 1;
for i = lowX:highX
    i;
    for j = lowY:highY
        for k = lowZ:highZ
%             position = (i-lowX)*(highX-lowX+1)^2 + (j-lowY)*(highY-lowY+1) + k-lowZ+1
            indices(position,:) = [i,j,k]; 
            position = position + 1;
        end 
    end
end
            
index = sub2ind(ct.cubeDim,indices(:,1),indices(:,2),indices(:,3));          
cst{1,4} = {index};   


lowX    = 140;
lowY    = 140;
lowZ    = 140;
highX   = 180;
highY   = 180;
highZ   = 180;

indices = zeros((highX-lowX+1)*(highY-lowY+1)*(highZ-lowZ+1),3);

position = 1;
for i = lowX:highX
    i;
    for j = lowY:highY
        for k = lowZ:highZ
%             position = (i-lowX)*(highX-lowX+1)^2 + (j-lowY)*(highY-lowY+1) + k-lowZ+1
            indices(position,:) = [i,j,k]; 
            position = position + 1;
        end 
    end
end
            
index = sub2ind(ct.cubeDim,indices(:,1),indices(:,2),indices(:,3));          
cst{2,4} = {index}; 


ct.cube{1} = ones(ct.cubeDim);
ct.cube{1}(1,1,1) = 0;

ct.cubeHU{1} = zeros(ct.cubeDim);
ct.cube{1}(1,1,1) = -1024;


ct.resolution.x = 0.5;
ct.resolution.y = 0.5;
ct.resolution.z = 0.5;

pln.propStf.isoCenter       = [ct.cubeDim(1) / 2 * ct.resolution.x, ...
                                ct.cubeDim(2) / 2 * ct.resolution.y, ...
                                ct.cubeDim(3) / 2 * ct.resolution.z];
 
 matRadGUI