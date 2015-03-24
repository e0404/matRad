function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mTissueClass,Interp,visBool)

[vRadDepthsSort, SortIndex] = sort(vRadDepths);
vRadDepthsSort= vRadDepthsSort./10;

vBeta = zeros(length(vRadDepths),1);
R0=sEnergy.range/10;


unsorted =1:length(vRadDepths);
newInd(SortIndex) =unsorted;


[~, index] = min(abs([Interp(:).energy]-sEnergy.energy));
NumTissueClass = size(Interp(1).alpha,2);
for i = 1:NumTissueClass
    mAlphaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).alpha(:,i), vRadDepthsSort,'pchip');
end

mAlpha=double(mAlphaSorted(newInd,:));
vAlpha=mAlpha(mTissueClass(:,2));


if min(Interp(index).beta(:)) == max(Interp(index).beta(:))
    vBeta(:) = 0.05;
else
    for i = 1:NumTissueClass
        mBetaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).beta(:,i), vRadDepthsSort,'pchip');
    end
    mBeta=double(mBetaSorted(newInd,:));
    vBeta=mBeta(mTissueClass(:,2));
end
%%
%mean = vAlphaSorted;
%std = vAlphaSorted.*0.2;
%vAlphaSorted = mean +std.*randn(numel(mean),1);
%figure,plot(vAlphaSortedNoise)

%% 
%vAlphaSorted = vAlphaSorted-vAlphaSorted.*0.25;

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







