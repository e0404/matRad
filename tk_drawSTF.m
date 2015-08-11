figure
hold on
for i=1:numel(stf(1).ray)
    plot(stf(1).ray(i).rayPos_bev(:,1),stf(1).ray(i).rayPos_bev(:,3),'xr')
end
axis([-106 106 -106 106])
axis ij