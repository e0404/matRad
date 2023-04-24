function [ThetaMatrix,ThetaVector] = matRad_getThetaMatrix(templateNormal,DistanceMatrix)
% getThetaMatrix gets (seed x dosepoint) matrix of relative polar angles
%
% call
%   [ThetaMatrix,ThetaVector] = matRad_getThetaMatrix(templateNormal,...
%       DistanceMatrix)
%   normally called within matRad_getBrachyDose
%   !!getDistanceMatrix needs to be called first!!
%
% input
%   DistanceMatrix:     [dosePoint x seedPoint] struct with fields 'x','y',
%                       'z' and total distance 'dist'
%   templateNormal:     normal vector of template (its assumed that this is
%                       the dir all seeds point to)
%
% output
%   angle matrix:       rows: index of dosepoint 
%                       columns: index of deedpoint
%                       entry: polar angles betreen seedpoints and  
%                       dosepoint in degrees
%   angle vector:       column vector of angle matrix entries
%
% comment:
%   The shape of the Theta matrix will be consistent with the shape of 
%   input fields. 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DistanceMatrix.dist(DistanceMatrix.dist == 0) = 1; %Avoid deviding by zero

ThetaMatrix = acosd((templateNormal(1)*DistanceMatrix.x + templateNormal(2)*DistanceMatrix.y + templateNormal(3)*DistanceMatrix.z)./DistanceMatrix.dist);  
if nargout == 2
    ThetaVector = reshape(ThetaMatrix,[],1);
end

end

