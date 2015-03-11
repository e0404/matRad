function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mT,Interp,visBool)



vAlpha = interp3(Interp.X,Interp.Y,Interp.Z,Interp.V,...
    (vRadDepths)./10,ones(size(vRadDepths,1),1)*sEnergy.energy,mT(:,2));
vAlpha(isnan(vAlpha))=0;
vBeta = 0.04;

%% plot interpolation

if visBool == 1
  [vRadDepthsSort, SortIndex] = sort(vRadDepths);
   vAlphaSort = vAlpha(SortIndex); 
   str = sprintf('interpolated alpha for %d MeV considering different tissues along z',round(sEnergy.energy));
   figure,plot(vRadDepthsSort./10,vAlphaSort,'LineWidth',2),title(str),...
      xlabel('radiological depth in cm'),ylabel('alpha'),grid minor;
end







