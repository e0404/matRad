function [ThetaMatrix,ThetaVector] = getThetaMatrix(templateNormal,DistanceMatrix)
%% seed dosepoint angle matrix
% !!getDistanceMatrix needs to be called first!!
% INPUT
% - DistanceMatrix [dosePoint x seedPoint] Struct with x,y,z and total distance fields
% - templateNorml: normal vector of template (its assumed that this is the dir all seeds point to)
% OUTPUT
% - angle matrix:
%       rows: index of dosepoint 
%       columns: index of deedpoint
%       entry: polar angles betreen seedpoints and dosepoint in degrees
% - angle vector:
%       column vector of angle matrix entries

ThetaMatrix = acosd((templateNormal(1)*DistanceMatrix.x + templateNormal(2)*DistanceMatrix.y + templateNormal(3)*DistanceMatrix.z)./DistanceMatrix.norm);  
ThetaVector = reshape(ThetaMatrix,[],1);

end

