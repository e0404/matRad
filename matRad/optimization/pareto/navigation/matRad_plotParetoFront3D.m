function matRad_plotParetoFront3D(retStruct)
    % matRad_plotParetoFront3D implements a function to visualize a 3d Pareto
    % surface
    % 
    %
    % References
    %   
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ps = retStruct.finds;
    if size(ps,2) ~= 3
        ME = MException('Matlab:WrongSize','Function only supports plotting of 3D surfaces!');
        throw(ME);
    end
    

    [k,facets] = matRad_ParetoSurfFromFacets(ps);
    figure
    trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor',[0.8 0.8 0.8])
    hold on 
    scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
            'MarkerFaceColor',[0 0 0])
    legend('approx. surface','Location','northwest')
    xlabel('f_1[a.u]')
    ylabel('f_2[a.u]')
    zlabel('f_3[a.u]')
    view(-45,5)
end