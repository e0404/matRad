dir = 'C:\Users\eric\Documents\PhD Project\XCAT\Data from XCAT';

origResolution.x = 0.85;
origResolution.y = 0.85;
origResolution.z = 1;

ct.motionVecX = cell(ct.numOfCtScen,1);
ct.motionVecY = cell(ct.numOfCtScen,1);
ct.motionVecZ = cell(ct.numOfCtScen,1);

for i = 2:ct.numOfCtScen
    fprintf('\nReading in motion vector file for frame %d of %d.\n',i,ct.numOfCtScen);
    
    fname = sprintf('VMATlungmotion_motion_vec_frame1_to_frame%d.txt',i);
    ct = matRad_parseMVF(ct,origResolution,fname,dir);
end

save LUNG_XCAT