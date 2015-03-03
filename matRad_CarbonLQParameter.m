function [vAlpha vBeta]= matRad_CarbonLQParameter(vRadDepths,vRadialDist_sq,sEnergy,mT,cst,Interp,visBool)


R0 = sEnergy.range/10;
E0 = sEnergy.energy;

vAlpha = Interp(vRadDepths./10,ones(size(vRadDepths,1),1)*E0,mT(:,2));

vBeta = ones(size(vRadDepths,1),1)*0.4;


vIsNan = isnan(vAlpha);
if(sum(vIsNan(:))>0)
    error('some alpha values couldnt be interpolated')
end












