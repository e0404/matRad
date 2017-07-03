% cubo = zeros([512,512,39]);
% cubo(ans) = L;
% cubo(idx_shift) = 6;
% cubo = zeros([160 160 160]);
% cubo(initIx) = 1;
% cubo(idx) = 10;

% cubo = ct.cube{1};
% cubo = rocubo;

% cubo = TracMatf;




% figure
hold off
scatter3(80, 80, 80,'xr')
hold on
surf([0,160;0,160],[0,0;0,0],[0,0;160,160],'FaceAlpha',0.1,'FaceColor','c')
surf([0,0;0,0],[0,160;0,160],[0,0;160,160],'FaceAlpha',0.1,'FaceColor','c')
surf([0,160;0,160],[160,160;160,160],[0,0;160,160],'FaceAlpha',0.1,'FaceColor','c')
surf([160,160;160,160],[0,160;0,160],[0,0;160,160],'FaceAlpha',0.1,'FaceColor','c')

% surf([-256,256;-256,256],[-256,-256;-256,-256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([-256,-256;-256,-256],[-256,256;-256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([-256,256;-256,256],[256,256;256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')
% surf([256,256;256,256],[-256,256;-256,256],[-195,-195;195,195],'FaceAlpha',0.1,'FaceColor','c')

% surf([stf.isoCenter(2)-512,stf.isoCenter(2);stf.isoCenter(2)-512,stf.isoCenter(2)],[stf.isoCenter(1)-512,stf.isoCenter(1)-512;stf.isoCenter(1)-512,stf.isoCenter(1)-512],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2)-512,stf.isoCenter(2)-512;stf.isoCenter(2)-512,stf.isoCenter(2)-512],[stf.isoCenter(1)-512,stf.isoCenter(1);stf.isoCenter(1)-512,stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2)-512,stf.isoCenter(2);stf.isoCenter(2)-512,stf.isoCenter(2)],[stf.isoCenter(1),stf.isoCenter(1);stf.isoCenter(1),stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')
% surf([stf.isoCenter(2),stf.isoCenter(2);stf.isoCenter(2),stf.isoCenter(2)],[stf.isoCenter(1)-512,stf.isoCenter(1);stf.isoCenter(1)-512,stf.isoCenter(1)],[stf.isoCenter(3)-390,stf.isoCenter(3)-390;stf.isoCenter(3),stf.isoCenter(3)],'FaceAlpha',0.1,'FaceColor','c')

% surf([0,512;0,512],[0,0;0,0],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([0,0;0,0],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([0,512;0,512],[512,512;512,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')
% surf([512,512;512,512],[0,512;0,512],[0,0;39,39],'FaceAlpha',0.1,'FaceColor','c')


view(3)

level = [1e-05,9,100];
face = [.3,.9,1];
colo = ['b','g','k'];
for i = 1:1
[f,v] = isosurface(cubo, level(i));
patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',face(i),'FaceColor',colo(i));
end


% iso = [80 80 80];
% iso = [stf.isoCenter(1),stf.isoCenter(2),stf.isoCenter(3)]./3;

% target2 = stf.ray(1).targetPoint./3 + iso;
% source2 = stf.sourcePoint./3 + iso;


line([target2(1) source2(1)],[target2(2) source2(2)],[target2(3) source2(3)])
% line([80 80],[-1000 1000],[80 80])
line([target3(1) source3(1)],[target3(2) source3(2)],[target3(3) source3(3)],'Color','r')
axis([-30 190 -30 190 0 160])

% scatter3(D(1:100:end,1),D(1:100:end,2),D(1:100:end,3),'xr')

% cubo2 = zeros([512,512,39]);
% cubo2(idx) = 1;
% 
% [f,v] = isosurface(cubo2, 0.9);
% patch('Faces',f,'Vertices',v,'EdgeAlpha',0.0001,'FaceAlpha',.8,'FaceColor','g');
