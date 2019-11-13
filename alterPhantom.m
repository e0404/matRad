matRad_rc

load BOXPHANTOM

 %% Input parameters

ct.resolution.x = 1;
ct.resolution.y = 1;
ct.resolution.z = 1;

lengthX = 100;
lengthY = 200;
lengthZ = 100;

targetCenter = [50, 50, 50];
targetWidth  = [5,  5, 5];

ct.cubeDim = [round(lengthY/ct.resolution.y), round(lengthX/ct.resolution.x), round(lengthZ/ct.resolution.z)];

lowX   = 1;
lowY   = 1;
lowZ   = 1;
highX  = ct.cubeDim(1);
highY  = ct.cubeDim(2);
highZ  = ct.cubeDim(3);


 %% Calculation of phantom

 
targetCenter = [targetCenter(2), targetCenter(1), targetCenter(3)];
targetWidth  = [targetWidth(2) , targetWidth(1) , targetWidth(3)];


ct.cube{1} = ones(ct.cubeDim);
ct.cube{1}(1,1,1) = 0;

ct.cubeHU{1} = zeros(ct.cubeDim);
ct.cube{1}(1,1,1) = -1024;



indices = zeros((highX-lowX+1)*(highY-lowY+1)*(highZ-lowZ+1),3);

position = 1;
for i = lowX:highX
    i
    for j = lowY:highY
        for k = lowZ:highZ
            indices(position,:) = [j,i,k]; 
            position = position + 1;
        end 
    end
end
            
index = sub2ind([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)],indices(:,1),indices(:,2),indices(:,3));          
cst{1,4} = {index};   


lowX    = round((targetCenter(1)-targetWidth(1)/2)/ct.resolution.x);
lowY    = round((targetCenter(2)-targetWidth(2)/2)/ct.resolution.y);
lowZ    = round((targetCenter(3)-targetWidth(3)/2)/ct.resolution.z);
highX   = round((targetCenter(1)+targetWidth(1)/2)/ct.resolution.x);
highY   = round((targetCenter(2)+targetWidth(2)/2)/ct.resolution.y);
highZ   = round((targetCenter(3)+targetWidth(3)/2)/ct.resolution.z);

indices = zeros((highX-lowX+1)*(highY-lowY+1)*(highZ-lowZ+1),3);

position = 1;
for i = lowX:highX
    i
    for j = lowY:highY
        for k = lowZ:highZ
            indices(position,:) = [i,j,k]; 
            position = position + 1;
        end 
    end
end
            
index = sub2ind(ct.cubeDim,indices(:,1),indices(:,2),indices(:,3));          
cst{2,4} = {index}; 

save('phantomTest.mat', 'ct', 'cst');


