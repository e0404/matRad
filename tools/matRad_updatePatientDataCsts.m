
load('BOXPHANTOM.mat');
save('BOXPHANTOM_old.mat');    
cst = matRad_convertOldCstToNewCstObjectives(cst);
save('BOXPHANTOM.mat');

load('HEAD_AND_NECK.mat');
save('HEAD_AND_NECK_old.mat');
cst = matRad_convertOldCstToNewCstObjectives(cst);
save('HEAD_AND_NECK.mat');

load('LIVER.mat');
save('LIVER_old.mat');
cst = matRad_convertOldCstToNewCstObjectives(cst);
save('LIVER.mat');

load('PROSTATE.mat');
save('PROSTATE_old.mat');
cst = matRad_convertOldCstToNewCstObjectives(cst);
save('PROSTATE.mat');

load('TG119.mat');
save('TG119_old.mat');
cst = matRad_convertOldCstToNewCstObjectives(cst);
save('TG119.mat');
    
    
    
    
 