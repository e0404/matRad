

function maskCube = matRad_CtDownsamplingMask(ct,res)
% getting new coarse physical coordinates
X = ct.x(1):res(1):ct.x(end);
Y = ct.y(1):res(2):ct.y(end);
Z = ct.z(1):res(3):ct.z(end);
% index of nearest points (X)
xx = (ct.x'*ones(1,numel(X))); 
XX = ((X)' * ones(1,numel(ct.x)))';
dif = xx-XX;
minval = min(abs(dif));
for i=1:size(dif,2), Xidx(i) = find(minval(i) == abs(dif(:,i)));end

% index of nearest points  (Y)
yy = (ct.y'*ones(1,numel(Y)));
YY = ((Y)' * ones(1,numel(ct.y)))';
dif = yy-YY;
minval = min(abs(dif));
for i=1:size(dif,2), Yidx(i) = find(minval(i) == abs(dif(:,i)));end

% index of nearest points  (Z)
zz = (ct.z'*ones(1,numel(Z)));
ZZ = ((Z)' * ones(1,numel(ct.z)))';
dif = zz-ZZ;
minval = min(abs(dif));
for i=1:size(dif,2), Zidx(i) = find(minval(i) == abs(dif(:,i)));end

maskCube = zeros(ct.cubeDim);

for i=1:numel(Xidx)
    for j=1:numel(Yidx)
        for k = 1: numel(Zidx)
            maskCube(Xidx(i), Yidx(j),Zidx(k)) = 1;
        end 
    end 
end 
end
