function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mT,Interp,visBool)

[vRadDepthsSort, SortIndex] = sort(vRadDepths);
vRadDepthsSort= vRadDepthsSort./10;


R0=sEnergy.range/10;

unsorted =1:length(vRadDepths);
newInd(SortIndex) =unsorted;


[~, index] = min(abs(Interp.vEnergies-sEnergy.energy));
vAlphaSorted = interp1(R0-Interp.mDepth(:,index), Interp.mAlpha(:,index), vRadDepthsSort,'pchip');

%%
%mean = vAlphaSorted;
%std = vAlphaSorted.*0.2;
%vAlphaSorted = mean +std.*randn(numel(mean),1);
%figure,plot(vAlphaSortedNoise)

%% 
%vAlphaSorted = vAlphaSorted-vAlphaSorted.*0.25;

vAlpha=double(vAlphaSorted(newInd));
vBeta = 0.05;

%
%str =sprintf('Range of this beam is %f',R0);
%figure,plot(vRadDepthsSort,vAlphaSorted),title(str);

% figure,subplot(121),plot(vRadDepthsSort,vAlphaSorted),title('own interpolation');
%        subplot(122),plot(vRadDepthsSort,vAlpha2),title('MTPS');
        



%% plot interpolation

if visBool == 1

   str = sprintf('interpolated alpha for %d MeV considering different tissues along z',round(sEnergy.energy));
   figure,plot(vRadDepthsSort./10,vAlphaSorted,'LineWidth',2),title(str),...
      xlabel('radiological depth in cm'),ylabel('alpha'),grid minor
  
  
end







