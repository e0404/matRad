function [k,facets] = matRad_ParetoSurfFromFacets(fVals)
% matRad_ParetoSurfFromFacets allows to calculate the facets of the pareto surface
% from its vertices. The facets are calculated by calculating the convex hull.
%
% call
%   [k,facets] = matRad_ParetoSurfFromFacets(fVals)
%
% input
%   fVals:          Matrix storing the objective function values of the pareto surface vertices (objective function values of the calculated points)
%
% output
%   k:              Matrix storing the vertices belonging to the facets of the convex hull
%   facets:         Matrix storing only the vertices that belong to the Pareto surface (The ones associated to facets with positive normals)
%
% References
%
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

    %calculate convex hull
    [k,vol] = convhulln(fVals);

    %%
    %calculate normals for each facet(normal vectors are perpendicular to their respective facet)
    %Done by solving P*n= 0 where P is a matrix with the vectors spanning the hyperplane and n the
    %normal vector to be calculated

    %initializing some objects that are returned by the function
    facets= zeros(size(k));

    j = 0;
    %% loop over all facets of convex hull
    for i = 1:size(k,1)
        %% Step2: Calculate upper bounds from convex hull
        %vertices of facets
        ps = fVals(k(i,:),:);
        %calculate vectors spanning hyperplane through refPoint
        refPoint = ps(1,:);
        f = ps-refPoint; % vectors spanning hyperplan through ps(1,:) 
        zw = f(2:end,:);  % remove reference Point
        %
        normal = null(zw); %solve P*n=0 P Matrix
    
        %% check orientation of calculated normal vector
        %get vertex not in facet and calculate orientation of facet
        idxs = (1:size(fVals,1));
        diffs = setdiff(idxs,k(i,:)); % find vertex not in facet
        
        %check if orientation and normal vector face in same direction
        orientationVector = fVals(diffs(1),:)-refPoint;
        orientation =(orientationVector*normal>0);    

        %flip orientation of normal vector if it goes in the opposite direction
        normal = normal*(2*orientation-1);
        normal = round(normal,4);
        
        %reject facet if normal has all negative components
        
        if any(round(normal,3)<0)
            continue
        end
        facets(i,:) = k(i,:);
    end
end