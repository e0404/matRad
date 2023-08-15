function matRad_plotParetoSurface(retStruct)
% matRad function that plots a color coded Pareto Surface. Colors are based
% on penalty values of the data points in 3 dimensions.
% Red means higher penalties of objective plotted along x-axis.
% Green means higher penalties of objective plotted along y-axis.
% Blue means higher penalties of objective plotted along z-axis.
%
% call
%   cBarHandle = matRad_plotColorbar(axesHandle,cMap,window,varargin)
%
% input
%   axesHandle  handle to axes the colorbar will be attached to
%   cMap        corresponding colormap
%   window      colormap window (corresponds to clim)
%   varargin    additional key-value pairs that will be forwarded to the
%               MATLAB colorbar(__) call
%
% output
%   cBarHandle  handle of the colorbar object
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
fig2 = figure;
switch size(fInd,2)
    case 2
        penGrid = [penGrid,zeros(size(penGrid,1),1)];
        scatter(fInd(:,1),fInd(:,2),[],penGrid,'filled')
        xlabel("x: " + VOIObj(1));
        ylabel("y: " + VOIObj(2));
        %set(gca,'Xscale','log')
        %set(gca,'Yscale','log') 
    case 3
        scatter3(fInd(:,1),fInd(:,2),fInd(:,3),[], penGrid,'filled')
        xlabel("x: " + VOIObj(1));
        ylabel("y: " + VOIObj(2));
        zlabel("z: " + VOIObj(3));
    otherwise
        warning(['Number of objectives for Pareto Analysis not suited for Plot!']);
end
%}
figure;
ps = retStruct.finds;
[k,facets] = matRad_ParetoSurfFromFacets(ps);

trisurf(facets(all(facets,2),:),ps(:,1),ps(:,2),ps(:,3),'FaceColor',[0.8 0.8 0.8])
hold on 
scatter3(ps(:,1),ps(:,2),ps(:,3),'MarkerEdgeColor','black',...
        'MarkerFaceColor',[0 0 0])
