function [DistanceMatrix,DistanceVector] = getDistanceMatrix(seedPoints,dosePoints)
%% get seed dosepoint distance matrix
% INPUT
% - seedPoints struct with fields x,y,z
% - dosePoints struct with fields x,y,z
% OUTPUT
% - distance matrix:
%       rows: index of dosepoint 
%       columns: index of deedpoint
%       entry: distance of seedpoints and dosepoint in cm
% - distance vector:
%       column vector of distance matrix entries

DistanceMatrix.x = dosePoints.x'*ones(1,length(seedPoints.x)) - ones(length(dosePoints.x),1)*seedPoints.x;
DistanceMatrix.y = dosePoints.y'*ones(1,length(seedPoints.y)) - ones(length(dosePoints.y),1)*seedPoints.y;
DistanceMatrix.z = dosePoints.z'*ones(1,length(seedPoints.z)) - ones(length(dosePoints.z),1)*seedPoints.z;
DistanceMatrix.norm = sqrt(DistanceMatrix.x.^2+DistanceMatrix.y.^2+DistanceMatrix.z.^2);
if nargout == 2
DistanceVector = reshape(DistanceMatrix.norm,[],1);
end

end

