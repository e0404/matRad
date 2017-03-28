function [gammaCube,myColormap,gammaPassRate] = matRad_gammaIndex(cube1,cube2,resolution,criteria,slice,interpoints,localglobal)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamma index calculation according to http://www.ncbi.nlm.nih.gov/pubmed/9608475
% 
% call
%    gammaCube = matRad_gammaIndex(cube1,cube2,resolution,slice)
%
% input
%   cube1:          dose cube as an M x N x O array
%   cube2:          dose cube as an M x N x O array
%   resolution:     resolution of the cubes [mm/voxel]
%   slice:          optional parameter for plotting the gamma pass rate on
%                   a specific slice
%   criteria:      optional [1x2] vector depicting the gamma criteria      
%   interpoints:   optional paramenter that represents the number of
%                   iterative interpolations on the input cubes (e.g. if 
%                   the input is k there will be 2^k-1 interpolation 
%                   points). The maximum suggested value is 3. 
%                   default value: 0
%                   Interpolation is made with cubic interpolation.
%   localglobal:    optional parameter that choose the kind of gamma want 
%                       to evaluate
%                   + 'global': default
%                   + 'local': 
%
% output 
%   gammaCube       relative distance between start and endpoint for the 
%
%   gamma index calculation represents a technique for 
%   the quantitative evaluation of dose distributions.
% References
%   [1]  http://www.ncbi.nlm.nih.gov/pubmed/9608475
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert
%
% This file is NOT part of the official matRad release. 
% This file has to be used only for internal purposes! 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if cubes consistent
if ~isequal(size(cube1),size(cube2))
   error('dose cubes must be the same size\n'); 
end

% set parameters for gamma index calculation
if exist('interpoints','var')
    n = interpoints;
else
    n = 0;
end

% set parameters for gamma index calculation
if exist('criteria','var')
    dist2AgreeMm     = criteria(1); % in [mm]
    relDoseThreshold = criteria(1); % in [%]
else
    dist2AgreeMm     = 1; % in [mm]
    relDoseThreshold = 1; % in [%]
end

% set global or local gamma evaluation
if ~exist('localglobal','var')
    localglobal = 'global';
end

%fprintf(['matRad: using gamme criteria ' num2str(dist2AgreeMm) '[mm], ' num2str(relDoseThreshold) '[%%]' ]);

% convert to absolute doses (use global max) and distance in voxels
absDoseThreshold = relDoseThreshold/100 * max(cube1(:));

% define search nighborhood and resolution
if round(n) ~= n
    error('n must be an integer value');
else
    resolution=resolution./(2.^n);
end

neighborhood = 1.5* dist2AgreeMm; % [mm]
searchX = ceil(neighborhood./resolution(1));
searchY = ceil(neighborhood./resolution(2));
searchZ = ceil(neighborhood./resolution(3));

% init cube
% cut all the zeros surrounding the interesting part
cut1 = any(any(cube1,2),3);
cut2 = any(any(cube1,1),3);
cut3 = any(any(cube1,1),2);

% avoids that "zero-slides" are deleted between two interest regions
k=find(cut1);
cut1(k(1):k(end))=1;
k=find(cut2);
cut2(k(1):k(end))=1;
k=find(cut2);
cut2(k(1):k(end))=1;

cubex1c=cube1;

cubex1c( ~cut1, : , : ) = [];  %rows
cubex1c( :, ~cut2 , : ) = [];  %columns
cubex1c( :, :, ~cut3) = [];  

% add some voxels with value zero around the cube
cubex1 = zeros([size(cubex1c,1)+2*searchX size(cubex1c,2)+2*searchY size(cubex1c,3)+2*searchZ]);
cubex1((1+searchX):(end-searchX), ...
            (1+searchY):(end-searchY), ...
            (1+searchZ):(end-searchZ)) = cubex1c;
        
cubex2c=cube2;

cubex2c( ~cut1, : , : ) = [];  %rows
cubex2c( :, ~cut2 , : ) = [];  %columns
cubex2c( :, :, ~cut3) = [];  

cubex2 = zeros([size(cubex2c,1)+2*searchX size(cubex2c,2)+2*searchY size(cubex2c,3)+2*searchZ]);
cubex2((1+searchX):(end-searchX), ...
            (1+searchY):(end-searchY), ...
            (1+searchZ):(end-searchZ)) = cubex2c;

% eventual interpolation
if n > 0
    %cubex1 = interp3(cubex1,n,method);
    cubex2 = interp3(cubex2,n,'cubic');
end

tmpCube = zeros(size(cubex1)); 
tmpCube((1+searchX):(end-searchX), ...
            (1+searchY):(end-searchY), ...
            (1+searchZ):(end-searchZ)) = inf;

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
                           - cubex2((1+((2^n)*searchX)+i) : 2^n : (end-((2^n)*searchX)+i), ...
                                   (1+((2^n)*searchY)+j) : 2^n : (end-((2^n)*searchY)+j), ...
                                   (1+((2^n)*searchZ)+k) : 2^n : (end-((2^n)*searchZ)+k));

            if exist('localglobal','var') && strcmp(localglobal,'local')
                absDoseThreshold = cubex1 .* relDoseThreshold/100;                
            end
                    
            tmpCube = tmpCube.^2 ./ absDoseThreshold.^2 + delta_sq;
            %gammaCubeSq(:) = min([gammaCubeSq(:) tmpCube(:)],[],2);
            gammaCubeSq(ix) = min(gammaCubeSq(ix),tmpCube(ix));
            
            
        end
    end
    
%    display '.';
    
end

% evaluate gamma cube and set to zero all the voxel that contain an
% infinite value
gammaCubeSqx = zeros(size(cube1));
gammaCubeSqx(cut1,cut2,cut3) = gammaCubeSq((1+searchX):(end-searchX), ...
            (1+searchY):(end-searchY), ...
            (1+searchZ):(end-searchZ));

gammaCube = sqrt(gammaCubeSqx);
infcube=isinf(gammaCube);
gammaCube(infcube)=0;

    
    myColormap = [  0.2081    0.1663    0.5292
                    0.2336    0.1932    0.5444
                    0.2592    0.2201    0.5596
                    0.2847    0.2470    0.5748
                    0.3103    0.2739    0.5899
                    0.3358    0.3008    0.6051
                    0.3614    0.3277    0.6203
                    0.3869    0.3546    0.6355
                    0.4125    0.3814    0.6507
                    0.4380    0.4083    0.6659
                    0.4636    0.4352    0.6811
                    0.4891    0.4621    0.6963
                    0.5146    0.4890    0.7114
                    0.5402    0.5159    0.7266
                    0.5657    0.5428    0.7418
                    0.5913    0.5697    0.7570
                    0.6168    0.5966    0.7722
                    0.6424    0.6235    0.7874
                    0.6679    0.6504    0.8026
                    0.6935    0.6773    0.8178
                    0.7190    0.7042    0.8329
                    0.7445    0.7311    0.8481
                    0.7701    0.7580    0.8633
                    0.7956    0.7849    0.8785
                    0.8212    0.8117    0.8937
                    0.8467    0.8386    0.9089
                    0.8723    0.8655    0.9241
                    0.8978    0.8924    0.9393
                    0.9234    0.9193    0.9544
                    0.9489    0.9462    0.9696
                    0.9745    0.9731    0.9848
                    1.0000    1.0000    1.0000
                    1.0000    1.0000         0
                    1.0000    0.9677         0
                    1.0000    0.9355         0
                    1.0000    0.9032         0
                    1.0000    0.8710         0
                    1.0000    0.8387         0
                    1.0000    0.8065         0
                    1.0000    0.7742         0
                    1.0000    0.7419         0
                    1.0000    0.7097         0
                    1.0000    0.6774         0
                    1.0000    0.6452         0
                    1.0000    0.6129         0
                    1.0000    0.5806         0
                    1.0000    0.5484         0
                    1.0000    0.5161         0
                    1.0000    0.4839         0
                    1.0000    0.4516         0
                    1.0000    0.4194         0
                    1.0000    0.3871         0
                    1.0000    0.3548         0
                    1.0000    0.3226         0
                    1.0000    0.2903         0
                    1.0000    0.2581         0
                    1.0000    0.2258         0
                    1.0000    0.1935         0
                    1.0000    0.1613         0
                    1.0000    0.1290         0
                    1.0000    0.0968         0
                    1.0000    0.0645         0
                    1.0000    0.0323         0
                    1.0000         0         0];
 
% need to modify the absDoseThreshold adding a 
if exist('localglobal','var') && strcmp(localglobal,'local')
    absDoseThresholdx = zeros(size(cube1));
    absDoseThresholdx(cut1,cut2,cut3) = absDoseThreshold((1+searchX):(end-searchX), ...
    (1+searchY):(end-searchY), (1+searchZ):(end-searchZ));
    absDoseThreshold = absDoseThresholdx;
end
        
        doseIx          = cube1 > absDoseThreshold;
        numOfPassGamma  = sum(gammaCube(doseIx) < 1);
        gammaPassRate   = 100 * numOfPassGamma / sum(doseIx(:));
        
if nargin > 4   
        figure
        set(gcf,'Color',[1 1 1]);
        imagesc(gammaCube(:,:,slice),[0 2])
        colormap(gca,myColormap);
        colorbar
        
        title({[num2str(gammaPassRate,5) '% of points > ' num2str(relDoseThreshold) '% pass gamma criterion (' num2str(relDoseThreshold) '% / ' num2str(dist2AgreeMm) 'mm)']; ['with ' num2str(2^n-1) ' interpolation points']});
end
