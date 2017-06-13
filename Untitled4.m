cubo = zeros([512,512,39]);
cubo(V(ix)) = 1;

figure
hold
surf([0,512;0,512],[0,0;0,0],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
surf([0,0;0,0],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
surf([0,512;0,512],[512,512;512,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
surf([512,512;512,512],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')

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
[f,v] = isosurface(cubo, 0.9);
patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',.1,'FaceColor','r');

iso = [stf.isoCenter(2),stf.isoCenter(1),stf.isoCenter(3)];

target2 = stf.ray(2).targetPoint_bev + iso;
source2 = stf.sourcePoint_bev + iso;

% line([target2(2) source2(2)],[target2(1) source2(1)],[target2(3)./10 source2(3)./10])
axis([-50 550 -50 550 0 39])

% scatter3(D(1:100:end,2),D(1:100:end,1),D(1:100:end,3)./10,'xr')

% cubo2 = zeros([512,512,39]);
% cubo2(idx) = 1;
% 
% [f,v] = isosurface(cubo2, 0.9);
% patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',.8,'FaceColor','g');
