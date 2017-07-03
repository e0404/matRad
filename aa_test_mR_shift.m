rad = 10;
[x,y,z] = sphere(50);
isosphere = [120 60 120];
initIx = sub2ind([160 160 160], round(reshape(y,[], 1).*rad+isosphere(2)),...
    round(reshape(x,[], 1).*rad+isosphere(1)), round(reshape(z,[], 1).*rad+isosphere(3)));

dim = [160 160 160];
res = [3 3 3];
iso = [240 240 240];
offs = [50 50 34];
Gang = 90;
Cang = 0;
source = 6850.*[sind(Gang).*cosd(Cang),-cosd(Gang).*cosd(Cang),-sind(Cang)] + offs;
target = -1.*source + offs;
Dx = 60;
Dz = 60;



[idx,idx_shift,target2, source2] = mR_shift(initIx,dim,source,target,iso, res, [Gang Cang], Dx, Dz);

target3 = target + iso./res;
source3 = source + iso./res;

aa_plot3ddose