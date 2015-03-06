function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,vRadialDist_sq,sEnergy,mT,Interp,mDesign,visBool)

% radial distance is not yet considered 
%% using scattered Interpolant
% e.g. 275 calls take 20 seconds for 3D Interpolation
% e.g. 743 calls take 20 seconds for 2D Interpolation
%  tic
%  vAlpha = Interp(sEnergy.range-vRadDepths,ones(size(vRadDepths,1),1)*sEnergy.energy,mT(:,2));
%  time1= toc;


%% using interp3 -> contains NaN values -> extrapolation is done beforehand in scatteredInterpolant
% in less than 20 seconds for all NumBixel calls
vAlpha = interp3(Interp.X,Interp.Y,Interp.Z,Interp.V,...
    sEnergy.range-vRadDepths,ones(size(vRadDepths,1),1)*sEnergy.energy,mT(:,2));

 
%% using griddata needs the same time as scattered interpolant but doesnt extrapolate -> NaN values

%% sorting the variables is not improving performance
% [vRadDepthsSort, SortIndex] = sort(sEnergy.range-vRadDepths);
% vT = mT(:,2);
% vAlphaSort = Interp(vRadDepthsSort,ones(size(vRadDepths,1),1)*sEnergy.energy,vT(SortIndex));
% newIdx(SortIndex) = 1:length(vRadDepths);
% vAlpha = vAlphaSort(newIdx);


vBeta = 0.04;
%vBeta = ones(size(vRadDepths,1),1)*0.4;


%  vIsNan = isnan(vAlpha);
%  if(sum(vIsNan(:))>0)
%      Disp('code goes here');
%      error('some alpha values couldnt be interpolated');
%  end
vAlpha(isnan(vAlpha))=0;

%% plot alpha
%  [vRadDepthsSort, SortIndex] = sort(vRadDepths);
%   vAlphaSort = vAlpha(SortIndex); 
%   str = sprintf('interpolated alpha for %f MeV considering different tissues along z',sEnergy.energy);
%   figure,plot(vRadDepthsSort,vAlphaSort),title(str),...
%      xlabel('radiological depth'),ylabel('alpha');









