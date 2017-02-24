function matRad_compareDoseCubes(cube1,cube2,resolution,cellName,slice,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comparison of 3D cubes
% 
% call
%   matRad_compareDoseCubes(cube1,cube2,resolution,filename)
%
% input
%   cube1:      first 3D array
%   cube2:      second 3D array
%   resolution: resolution
%   filename:   optional: name of ps file for output
%
% output
%
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is not part of matRad.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
   slice = 100; 
end
    
defaultLineWidth = 1.5;

% necessary for sobp evaluation
%cube2 = permute(cube2,[3 1 2]);

% downsample cubes
% cube1 = matRad_downsampleImageStack(cube1,0.5);
% cube2 = matRad_downsampleImageStack(cube2,0.5);

if ~isequal(size(cube1),size(cube2))
    error('cube dimensions inconsistent');
end

%% integral dose
relIntDoseDif = (1-sum(cube1(:))/sum(cube2(:)))*100;

fprintf(['Relative difference in integral dose: ' num2str(relIntDoseDif) '%%\n']);

figure
set(gcf,'Color',[1 1 1]);
%% transversal slices
subplot(3,2,1)
imagesc(cube1(:,:,slice))
colorbar
title([ cellName{1,1} ': slice ' num2str(slice)])

subplot(3,2,2)
imagesc(cube2(:,:,slice))
colorbar
title([ cellName{1,2} ': slice ' num2str(slice)])

norm = max(max(cube1(:,:,slice)));

ax3 = subplot(3,2,3);
imagesc(100*(cube1(:,:,slice)-cube2(:,:,slice))./norm)
myMap = matRad_getCostumColorbarDiff(cube1,cube2,slice,3);
colormap(ax3,myMap); colorbar;
title(['rel diff [%] ' cellName{1,1} ' - ' cellName{1,2} ' : slice ' num2str(slice)])

%% idds

x = resolution.x*[1/2:1:size(cube1,2)-1/2];

cube1IDD = sum(sum(cube1,3),1);
cube2IDD = sum(sum(cube2,3),1);

subplot(3,2,4)
hold on
plot(x,cube1IDD,'r','LineWidth',defaultLineWidth)
plot(x,cube2IDD,'b--','LineWidth',defaultLineWidth)
title(['IDD - rel diff int. dose ' num2str(relIntDoseDif)])
box on
grid on
legend(cellName,'location','northwest')

%% profile along depth
AxialSliceCube1 = squeeze(cube1(:,:,slice));
AxialSliceCube2 = squeeze(cube2(:,:,slice));

[A,CentralIx]   = max(sum(AxialSliceCube1,2));

cube1CentralRayProf = AxialSliceCube1(CentralIx,:);
cube2CentralRayProf = AxialSliceCube2(CentralIx,:);

subplot(3,2,5)
hold on
plot(x,cube1CentralRayProf,'r','LineWidth',defaultLineWidth)
plot(x,cube2CentralRayProf,'b--','LineWidth',defaultLineWidth)
box on
grid on
title('depth profile')
legend(cellName,'location','northwest')

%% profile along lateral direction entrance

y = resolution.y*[1/2:1:size(cube1,1)-1/2];

LatIndex = 510;
cube1LatProfileEnt = AxialSliceCube1(:,LatIndex);
cube2LatProfileEnt = AxialSliceCube2(:,LatIndex);

NormFac = 1;%max(cube1LatProfileEnt);
cube1LatProfileEnt = cube1LatProfileEnt/NormFac;
cube2LatProfileEnt = cube2LatProfileEnt/NormFac;
 
subplot(3,2,6)
hold on
plot(y,cube1LatProfileEnt,'r','LineWidth',defaultLineWidth)
plot(y,cube2LatProfileEnt,'b--','LineWidth',defaultLineWidth)
box on
grid on
title('lateral entrance profile')
legend(cellName,'location','northwest')

if nargin > 5
    annotation('textbox', [0 0.9 1 0.1],'String', filename, ...
    'EdgeColor', 'none','HorizontalAlignment', 'center','FontSize',16)
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpsc','-r300','-append',filename)
    %print('-dpdf','-r300',filename)    
end
    
%% calculate gamma pass rate
gammaCube = matRad_gammaIndex(cube1,cube2,[resolution.x resolution.y resolution.z],slice);

%% relative differences
figure,set(gcf,'Color',[1 1 1]);
subplot(3,2,1)
plot(x,100*(cube1IDD-cube2IDD)/max(cube2IDD),'r')
box on
grid on
title('rel diff [%] idd')

subplot(3,2,2)
plot(x,100*(cube1CentralRayProf-cube2CentralRayProf)/max(cube2CentralRayProf),'r')
box on
grid on
title('rel diff [%] central profile')

subplot(3,2,3)
plot(y,100*(cube1LatProfileEnt-cube2LatProfileEnt)/max(cube2LatProfileEnt),'r')
box on
grid on
title('rel diff [%] lateral profile @ entrance')

depthIx = 230;

cube1LatProfileIx = squeeze(cube1(:,depthIx,slice));
cube2LatProfileIx = squeeze(cube2(:,depthIx,slice));

subplot(3,2,4)
plot(y,100*(cube1LatProfileIx-cube2LatProfileIx)/max(cube2LatProfileIx),'r')
box on
grid on
title('rel diff [%] lateral profile before peak')

if nargin > 4
    annotation('textbox', [0 0.9 1 0.1],'String', filename, ...
    'EdgeColor', 'none','HorizontalAlignment', 'center')
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    print('-dpsc','-r300','-append',filename)
    
end




end