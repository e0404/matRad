function [vAlpha, vBeta]= matRad_CarbonLQParameter(vRadDepths,sEnergy,mTissueClass,Interp,visBool)

[vRadDepthsSort, SortIndex] = sort(vRadDepths);
vRadDepthsSort= vRadDepthsSort./10;

vBeta = zeros(length(vRadDepths),1);
R0=sEnergy.range/10;


unsorted =1:length(vRadDepths);
newInd(SortIndex) =unsorted;

if sEnergy.energy>230 && sEnergy.energy<250  
    str = 'asdf';
end


[~, index] = min(abs([Interp(:).energy]-sEnergy.energy));
NumTissueClass = size(Interp(1).alpha,2);
for i = 1:NumTissueClass
    mAlphaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).alpha(:,i), vRadDepthsSort,'pchip');
end

mAlpha=double(mAlphaSorted(newInd,:));


for i = 1:NumTissueClass
    mBetaSorted(:,i) = interp1(R0-Interp(index).res_range(:,i), Interp(index).beta(:,i), vRadDepthsSort,'pchip');
end

mBeta=double(mBetaSorted(newInd,:));

if size(mAlpha,2)>1
    vAlpha=zeros(length(mAlpha),1);
    vBeta=zeros(length(mBeta),1);
    for i = 1:length(mAlpha)
        vAlpha(i)=mAlpha(i,mTissueClass(i,2));
        vBeta(i)=mBeta(i,mTissueClass(i,2));
    end
else
    vAlpha=mAlpha;
    vBeta=mBeta;
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







