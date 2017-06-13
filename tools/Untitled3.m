tic;
figure
profile on
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],...
    [3 3],round(pln.isoCenter(1,3)/ct.resolution.z),'linear',3,'global');

profile off
profile viewer

toc

tic;

matRad_gammaIndex_old(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z));

toc

n=0;
% eventual interpolation
if exist('method','var')
    if ~strcmp(method,'standard')
        n=str2double(method(end)); % improvement in resolution
        dim = size(cubex2);
        [Xq,Yq,Zq] = meshgrid(1:(n+1)*dim(2),1:(n+1)*dim(1),1:(n+1)*dim(3));
        cubex2 = interp3(cubex2,Xq,Yq,Zq);
        %cubex1 = interp3(cubex1,str2double(n),method(1:end-1));
        %cubex2 = interp3(cubex2,log2(n+1),method(1:end-1)); 
        resolution=resolution./(n+1);
    end
end

tmpCube = inf(size(cubex1)); 
ix = cubex1 > 0;

gammaCubeSq = zeros(size(cubex1));
gammaCubeSq((1+searchX):(end-searchX), ...
            (1+searchY):(end-searchY), ...
            (1+searchZ):(end-searchZ)) = inf;


% search for min
for i = -searchX:searchX
    for j = -searchY:searchY
        for k = -searchZ:searchZ
            
            delta_sq = ((i*resolution(1))^2 + ...
                        (j*resolution(2))^2 + ...
                        (k*resolution(3))^2) / dist2AgreeMm^2;                 
            
            tmpCube((1+searchX):(end-searchX), ...
                    (1+searchY):(end-searchY), ...
                    (1+searchZ):(end-searchZ)) = ...
                             cubex1((1+searchX):(end-searchX), ...
                                   (1+searchY):(end-searchY), ...
                                   (1+searchZ):(end-searchZ)) ...
                           - cubex2(((1+searchX)+(n*searchX)+i) : n+1 : (end-searchX-(n*searchX)+i), ...
                                   ((1+searchY)+(n*searchY)+j) : n+1 : (end-searchY-(n*searchY)+j), ...
                                   ((1+searchZ)+(n*searchZ)+k) : n+1 : (end-searchZ-(n*searchZ)+k));
