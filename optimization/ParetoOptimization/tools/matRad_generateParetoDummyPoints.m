function fValsMod = matRad_generateParetoDummyPoints(fVals,U) 
    % matRad Function that generates Dummy Points for a set of points
    %
    % input
    %   fVals:      The so-far determined pareto optimal points 
    %   U:          Array storing highest objective values (should be an
    %               array of ones if normalized). Should be calculated after anchorpoint calculation and then kept constant to avoid issues from optimization complications
    %
    % output
    %   fValsMod:   Modified pareto set to remove mixed normal facets (see
    %               Rennen)
    %
    % References
    %   - https://dx.doi.org/10.2139/ssrn.1427721 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2023 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
theta = 0.1;
m = size(fVals,2);
n = size(fVals,1);

U = U*m+theta; %modifying upper bound 

fValsMod = zeros(size(fVals,1)*size(fVals,2),size(fVals,2)); %initialize

for i = 1:n %loop over and generate a dummy point for each dimension
    temp = repmat(fVals(i,:),m,1);
    for j = 1:numel(U)
        temp(j,j) = U(j);
    end
    fValsMod((i-1)*m+1:i*m,:) = temp;
end

fValsMod = [fVals;fValsMod];