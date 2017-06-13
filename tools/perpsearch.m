figure
hold
A = sourcePoint_bev;
B = targetPoint_bev;
line([sourcePoint_bev(2),targetPoint_bev(2)],[sourcePoint_bev(1),targetPoint_bev(1)],[sourcePoint_bev(3),targetPoint_bev(3)])
P = [33, 90, 12 ; 40, 95, 30];
scatter3(P(1,2), P(1,1), P(1,3))
scatter3(P(2,2), P(2,1), P(2,3))

Avec = repmat(A,2,1);
Bvec = repmat(B,2,1);

perpDist =  sqrt(sum(cross(P-Avec,P-Bvec).^2,2))./sqrt(sum((Bvec-Avec).^2,2));

otherDist = sqrt( sqrt(sum((P-Bvec).^2,2)).^2 - perpDist.^2 );

theta = asin(sqrt(sum((A(3)-B(3)).^2,2))/sqrt(sum((A-B).^2,2)));
phi = acos(sqrt(sum((A(1)-B(1)).^2,2))/sqrt(sum((A(1:2)-B(1:2)).^2,2)));

D = Bvec + [ otherDist*cos(theta)*cos(phi) -otherDist*cos(theta)*sin(phi) otherDist*sin(theta)];
scatter3(D(:,2),D(:,1),D(:,3))

