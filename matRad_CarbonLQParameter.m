function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mTissueClass,Interp,visBool)

% sort values for visualization and faster lookup
[vRadDepthsSort, SortIndex] = sort(vRadDepths);
vRadDepthsSort= vRadDepthsSort./10;

R0=sEnergy.range/10;
% remeber old indicies to desort them
unsorted =1:length(vRadDepths);
newInd(SortIndex) =unsorted;

% find corresponding index
[~, index] = min(abs([Interp(:).energy]-sEnergy.energy));
NumTissueClass = size(Interp(1).alpha,2);

for i = 1:NumTissueClass
    mAlphaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).alpha(:,i), vRadDepthsSort,'pchip');
end
% unsort values
mAlpha=double(mAlphaSorted(newInd,:));


for i = 1:NumTissueClass
    mBetaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).beta(:,i), vRadDepthsSort,'pchip');
end
% unsort values
mBeta=double(mBetaSorted(newInd,:));

% considering tissue specific alpha and beta values
linearInd = sub2ind(size(mAlpha),1:1:size(mAlpha,1),mTissueClass(:,2)');
vAlpha=mAlpha(linearInd)';
vBeta=mBeta(linearInd)';


%% plot interpolation
if visBool == 1

   str = sprintf('interpolated alpha for %d MeV considering different tissues along z',round(sEnergy.energy));
   figure,plot(vRadDepthsSort./10,vAlphaSorted,'LineWidth',2),title(str),...
      xlabel('radiological depth in cm'),ylabel('alpha'),grid minor
      waitforbuttonpress
  
end







