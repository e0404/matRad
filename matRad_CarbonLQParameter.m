function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mT,Interp,visBool)

[vRadDepthsSort, SortIndex] = sort(vRadDepths);
vRadDepthsSort= vRadDepthsSort./10;


R0=sEnergy.range/10;

unsorted =1:length(vRadDepths);
newInd(SortIndex) =unsorted;



%% rough interpolation
% tDepth = zeros(81,4);
% tDepth(:,1)=R0-Interp.tDepth(:,1);
% tDepth(:,2)=R0-Interp.tDepth(:,2);
% tDepth(:,3)=R0-Interp.tDepth(:,3);
% tDepth(:,4)=R0-Interp.tDepth(:,4);
% [c index] = min(abs(Interp.tTEnergies-sEnergy.energy));
% closestValues = Interp.tTEnergies(index);
% vAlphaSorted = interp1(tDepth(:,index), Interp.tTAlpha(:,index), vRadDepthsSort);
%% fine interpolation

  [c index] = min(abs(Interp.tTEnergies-sEnergy.energy));
  closestValues = Interp.tTEnergies(index);
  vAlphaSorted = interp1(R0-Interp.tDepth(:,index), Interp.tTAlpha(:,index), vRadDepthsSort,'pchip');

 %%
 
%mean = vAlphaSorted;
%std = vAlphaSorted.*0.2;

%vAlphaSorted = mean +std.*randn(numel(mean),1);

%figure,plot(vAlphaSortedNoise)

% add white Gaussian noise to signal

%vAlphaSorted = vAlphaSorted-vAlphaSorted.*0.25;

% for IX = 1 : numel(vRadDepthsSort)
%    
%         dummyAlpha = zeros(numel(Interp.tTEnergies),1);
%         
%         for JX = 1 : numel(Interp.tTEnergies)
%             dummyAlpha(JX) = interp1(R0-Interp.tDepth(:,JX), Interp.tTAlpha(:,JX), vRadDepthsSort(IX));
%         end
%         
%         vAlpha2(IX) = interp1(Interp.tTEnergies, dummyAlpha, sEnergy.energy);
%     
% end


vAlpha=vAlphaSorted(newInd);
%
%str =sprintf('Range of this beam is %f',R0);
%figure,plot(vRadDepthsSort,vAlphaSorted),title(str);

% figure,subplot(121),plot(vRadDepthsSort,vAlphaSorted),title('own interpolation');
%         subplot(122),plot(vRadDepthsSort,vAlpha2),title('MTPS');
        


% vAlpha = interp3(Interp.X,Interp.Y,Interp.Z,Interp.V,...
%     (vRadDepths)./10,ones(size(vRadDepths,1),1)*sEnergy.energy,mT(:,2));


% if ~all(vAlpha(:)) || sum(isnan(vAlpha))>0
%    str = 'some alpha values couldnt be interpolated'; 
% end

% 
% sigma = 0.1; 
% vAlphaNoise = vAlpha + sigma*randn(size(vAlpha));



% vAlpha(isnan(vAlpha))=0.54;
vBeta = 0.05;

%% plot interpolation

if visBool == 1
  [vRadDepthsSort, SortIndex] = sort(vRadDepths);
   vAlphaSort = vAlpha(SortIndex); 
   str = sprintf('interpolated alpha for %d MeV considering different tissues along z',round(sEnergy.energy));
   figure,plot(vRadDepthsSort./10,vAlphaSort,'LineWidth',2),title(str),...
      xlabel('radiological depth in cm'),ylabel('alpha'),grid minor;
  
  
%     [vRadDepthsSort, SortIndex] = sort(vRadDepths);
%    vAlphaSort = vAlpha(SortIndex); 
%    vAlphaSortNoise = vAlphaNoise(SortIndex); 
%    str = sprintf('interpolated alpha for %d MeV considering different tissues along z',round(sEnergy.energy));
%    figure,subplot(121),plot(vRadDepthsSort./10,vAlphaSort,'LineWidth',2),title(str),...
%       xlabel('radiological depth in cm'),ylabel('alpha'),grid minor;
%   subplot(122),plot(vRadDepthsSort./10,vAlphaSortNoise,'LineWidth',2),title(str),...
%       xlabel('radiological depth in cm'),ylabel('alpha'),grid minor;
  
  
  
end







