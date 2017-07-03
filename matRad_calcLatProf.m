function [latProf] = matRad_calcLatProf(baseData,rad_distancesSq,sig,dim,radDepths,ix)



% Here there is an approximation on radDepths. We use it the one of a
% test-ray. Will improve it in future.
X = matRad_interp1(baseData.depths,[baseData.sigma1 baseData.weight baseData.sigma2], radDepths);

sigmaSq_Narr = X(:,1).^2 + sig^2;
sigmaSq_Bro  = X(:,3).^2 + sig^2;

% calculate lateral profile
L_Narr =  exp( -rad_distancesSq ./ (2*sigmaSq_Narr))./(2*pi*sigmaSq_Narr);
L_Bro  =  exp( -rad_distancesSq ./ (2*sigmaSq_Bro ))./(2*pi*sigmaSq_Bro );
L = baseData.LatCutOff.CompFac * ((1-(X(:,2))).*L_Narr) + (X(:,2).*L_Bro);

% Need to be done better, this is only a test
bins = [0,20,48,56];
peakvec = find(baseData.Z >= baseData.Z(bins(4)+1));
bins(5) = peakvec(end);
tailvec = find(baseData.Z < baseData.Z(bins(1)+1));

% esce una merda dal calcolo degli indici. rivedere
profileLat = zeros(dim);
profileLat(ix) = L;


cut1 = any(any(profileLat,2),3);
cut2 = any(any(profileLat,1),3);
cut3 = any(any(profileLat,1),2);

% avoids that "zero-slides" are deleted between two interest regions
k = find(cut1); cut1(k(1):k(end)) = 1;
k = find(cut2); cut2(k(1):k(end)) = 1;
k = find(cut2); cut2(k(1):k(end)) = 1;

% cut cubes
profileLat( ~cut1, :, :) = []; % rows
profileLat( :, ~cut2, :) = []; % columns
profileLat( :, :, ~cut3) = []; % slices  

xdist = size(profileLat,1);
xdistZ = size(baseData.Z,1);

for i = 1:4
latProf.section(i) = [baseData.Z(bins(i)+1) baseData.Z(bins(i+1))];
latProf.Mask(i) = profileLat(floor(xdist*bins(i)/xdistZ)+1 : floor(xdist*bins(i+1)/xdistZ)+1, :,:);
latProf.Mask(i) = sum(latProf.Mask(1),1)./size(latProf.Mask(1),1);
end
latProf.Mask(5) = profileLat(end,:,:);
latProf.section(5) = [0 tailvec(1)];

disp(':)')


