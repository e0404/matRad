function matRad_plotPS(refObj)
    %{
    figure
    scatter3(finds(:,1),finds(:,2),finds(:,3))
    hold on
    scatter3(fNew(:,1),fNew(:,2),fNew(:,3),'filled')
    %}
    ps = refObj.fIndsAll;
    psRed = refObj.fIndsRed;
    fNew = refObj.fRef;
    [k,facets] = matRad_ParetoSurfFromFacets(ps);
    figure
    trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor',[0.8 0.8 0.8])
    hold on 
    scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
            'MarkerFaceColor',[0 0 0])
    scatter3(psRed(:,1),psRed(:,2),psRed(:,3),'MarkerEdgeColor','black',...
            'MarkerFaceColor',[1 1 0])
    scatter3(fNew(:,1),fNew(:,2),fNew(:,3),'filled','MarkerFaceColor','red')
    legend('approx. surface','Blocked areas','Available areas','Current point','Location','northwest')
    xlabel('f_1[a.u]')
    ylabel('f_2[a.u]')
    zlabel('f_3[a.u]')
    view(-45,5)
end