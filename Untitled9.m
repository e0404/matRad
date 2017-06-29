dosedose = full(doseTmpContainer{1});
dosedose = reshape(dosedose,dij.dimensions);
hold off
hold on
% surf([0,512;0,512],[0,0;0,0],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([0,0;0,0],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([0,512;0,512],[512,512;512,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([512,512;512,512],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')

% surf([-256,256;-256,256],[-256,-256;-256,-256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([-256,-256;-256,-256],[-256,256;-256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([-256,256;-256,256],[256,256;256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([256,256;256,256],[-256,256;-256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')

% surf([stf.isoCenter(2)-512,stf.isoCenter(2);stf.isoCenter(2)-512,stf.isoCenter(2)],[stf.isoCenter(1)-512,stf.isoCenter(1)-512;stf.isoCenter(1)-512,stf.isoCenter(1)-512],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2)-512,stf.isoCenter(2)-512;stf.isoCenter(2)-512,stf.isoCenter(2)-512],[stf.isoCenter(1)-512,stf.isoCenter(1);stf.isoCenter(1)-512,stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2)-512,stf.isoCenter(2);stf.isoCenter(2)-512,stf.isoCenter(2)],[stf.isoCenter(1),stf.isoCenter(1);stf.isoCenter(1),stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2),stf.isoCenter(2);stf.isoCenter(2),stf.isoCenter(2)],[stf.isoCenter(1)-512,stf.isoCenter(1);stf.isoCenter(1)-512,stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')

view(3)
%cubo = interp3(
maxi = max(max(max(dosedose)));
isovec = [0.001*maxi, 0.2*maxi, 0.5*maxi, 0.8*maxi, 0.95*maxi];
facevec = [0.1, 0.1, 0.2, 0.4, 1];
colorvec = ['c', 'b', 'g', 'y', 'r'];

for i=1:5
    [f,v] = isosurface(dosedose, isovec(i));
    patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',facevec(i),'FaceColor',colorvec(i));
end
