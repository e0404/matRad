figure
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z));


figure
subplot(2,2,1)
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z),'linear2');

subplot(2,2,2)
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z),'cubic2');

subplot(2,2,3)
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z),'linear4');

subplot(2,2,4)
matRad_gammaIndex(dose_5mm,dose_3mm,[ct.resolution.x ct.resolution.y ct.resolution.z],[3 3],round(pln.isoCenter(1,3)/ct.resolution.z),'cubic4');
