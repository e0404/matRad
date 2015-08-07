function tk_visSingleShape(tempShape,Ix_l, Ix_r,titleString,xMin,xMax)

mode = 'physical'; % 'index';

if ~exist('titleString','var')
    titleString = 'fluence Mask';
end

% 1. visualize binary mask
% figure
set(gcf,'position',[10 500 600 400])
imagesc(tempShape)
title(titleString,'Fontsize',14)

% 2. visualize leaf positions
if strcmp(mode,'index')
    xMax = size(tempShape,2);
    figure
    set(gcf,'position',[650 500 600 400])
    hold on
    axis ij
    axis([0.5 xMax+0.5 0.5 size(tempShape,1)+0.5])
    set(gca,'Color',[0.5 0 0]);
    for i=1:numel(Ix_l)
        fill([0.5 Ix_l(i)+0.5 Ix_l(i)+0.5 0.5],[i-1/2 i-1/2 i+1/2 i+1/2],'b')
        fill([Ix_r(i)+0.5 xMax+0.5 xMax+0.5 Ix_r(i)+0.5],[i-1/2 i-1/2 i+1/2 i+1/2],'b')    
    end
    axis tight
end

if strcmp(mode,'physical')
    if ~exist('xMax') || ~exist('xMin')
        xMin = min(Ix_l);
        xMax = max(Ix_r);
    end
    figure
    set(gcf,'position',[650 500 600 400])
    hold on
    axis ij
    axis([xMin xMax 0.5 size(tempShape,1)+0.5])
    set(gca,'Color',[0.5 0 0]);
    for i=1:numel(Ix_l)
        fill([xMin Ix_l(i) Ix_l(i) xMin],[i-1/2 i-1/2 i+1/2 i+1/2],'b')
        fill([Ix_r(i) xMax xMax Ix_r(i)],[i-1/2 i-1/2 i+1/2 i+1/2],'b')    
    end
    axis tight
end




end