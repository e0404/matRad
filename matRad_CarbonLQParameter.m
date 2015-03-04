function [vAlpha vBeta]= matRad_CarbonLQParameter(vRadDepths,vRadialDist_sq,sEnergy,mT,Interp,visBool)


R0 = sEnergy.range/10;

load('GSI_Chardoma_Carbon_BioData.mat');
    % get number measured points - determined through experiments 
    sVectorLength = 0;
    for i=1:numel(stBioData)
        for j=1:size(stBioData{1,i},2)       
            sVectorLength =sVectorLength + size(stBioData{1,i}(j).Depths,1);
        end    
    end

    mDesign = zeros(sVectorLength,4);
    idx = 1;

    for i=1:numel(stBioData)
        for j=1:size(stBioData{1,i},2)  
            tmpLength = size(stBioData{1,i}(j).Depths,1);
            %mDesign(idx:idx+tmpLength-1,1)=stBioData{1,i}(j).Depths;
            % R0 ist Zieltiefe des spots
            mDesign(idx:idx+tmpLength-1,1)=R0 - stBioData{1,i}(j).Depths; 
            mDesign(idx:idx+tmpLength-1,2)=stBioData{1,i}(j).Energy;
            mDesign(idx:idx+tmpLength-1,3)=i;
            mDesign(idx:idx+tmpLength-1,4)=stBioData{1,i}(j).Alpha;
            idx = idx+tmpLength;
        end    
    end
    
    Interp = scatteredInterpolant(mDesign(:,1:3),mDesign(:,4));




vAlpha = Interp(vRadDepths./10,ones(size(vRadDepths,1),1)*sEnergy.energy,mT(:,2));
vBeta = ones(size(vRadDepths,1),1)*0.4;


vIsNan = isnan(vAlpha);
if(sum(vIsNan(:))>0)
    error('some alpha values couldnt be interpolated')
end


  [vRadDepthsSort, SortIndex] = sort(vRadDepths);
  vAlphaSort = vAlpha(SortIndex);
  
  %figure,plot(vRadDepthsSort,vAlphaSort);









