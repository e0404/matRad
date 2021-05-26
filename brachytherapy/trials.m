% simple dose influence matrix with 1/r^2:
dij.physicalDose = zeros(length(dij.xPos),length(stf.seeds.xPos));
for i = 1:length(dij.xPos)
    for j = 1:length(stf.seeds.xPos)
        r = sqrt(...
            (stf.seeds.xPos(j)-dij.xPos(i))^2+...
            (stf.seeds.yPos(j)-dij.yPos(i))^2+...
            (stf.seeds.zPos(j)-dij.zPos(i))^2);
        
        dij.physicalDose(i,j) = 1/r^2;
    end
end

x=1:0.1:4
y=1:0.1:4
[X,Y]=meshgrid(x,y)
Z=sin(X).^2+cos(Y).^2
surf(X,Y,Z)
view(2)