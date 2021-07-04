function [ThetaMatrix,ThetaVector] = getThetaMatrix(templateNormal,DistanceMatrix)
% getThetaMatrix gets (seed x dosepoint) matrix of relative polar angles
%
% call
%   [ThetaMatrix,ThetaVector] = getThetaMatrix(templateNormal,DistanceMatrix)
%   normally called within matRad_getBrachyDose
%   !!getDistanceMatrix needs to be called first!!
%
% input
% - DistanceMatrix [dosePoint x seedPoint] Struct with x,y,z and total distance fields
% - templateNorml: normal vector of template (its assumed that this is the dir all seeds point to)
%
% output
% - angle matrix:
%       rows: index of dosepoint 
%       columns: index of deedpoint
%       entry: polar angles betreen seedpoints and dosepoint in degrees
% - angle vector:
%       column vector of angle matrix entries

ThetaMatrix = acosd((templateNormal(1)*DistanceMatrix.x + templateNormal(2)*DistanceMatrix.y + templateNormal(3)*DistanceMatrix.z)./DistanceMatrix.dist);  
ThetaVector = reshape(ThetaMatrix,[],1);

end

