x = DoseObjectives.matRad_SquaredDeviation;
y = DoseObjectives.matRad_SquaredOverdosing;

load('BOXPHANTOM.mat');
save('BOXPHANTOM_copy.mat');    
    cst{1,6} = {y};
        cst{1,6}{1,1}.parameters = {5};
        cst{1,6}{1,1}.penalty = 100;
    cst{2,6} = {x};
        cst{2,6}{1,1}.parameters = {60};
        cst{2,6}{1,1}.penalty = 800;
    save('BOXPHANTOM.mat');

load('HEAD_AND_NECK.mat');
save('HEAD_AND_NECK_copy.mat');
    cst{13,6} = {y};
        cst{13,6}{1,1}.parameters = {25};
        cst{13,6}{1,1}.penalty = 100;
    cst{14,6} = {y};
        cst{14,6}{1,1}.parameters = {25};
        cst{14,6}{1,1}.penalty = 100;
    cst{15,6} = {x};
        cst{15,6}{1,1}.parameters = {63};
        cst{15,6}{1,1}.penalty = 1000;
    cst{16,6} = {x};
        cst{16,6}{1,1}.parameters = {70};
        cst{16,6}{1,1}.penalty = 1000;
    cst{17,6} = {y};
        cst{17,6}{1,1}.parameters = {30};
        cst{17,6}{1,1}.penalty = 800;
save('HEAD_AND_NECK.mat');

load('LIVER.mat');
save('LIVER_copy.mat');
    cst{15,6} = {y};
        cst{15,6}{1,1}.parameters = {25};
        cst{15,6}{1,1}.penalty = 300;
    cst{16,6} = {x};
        cst{16,6}{1,1}.parameters = {45};
        cst{16,6}{1,1}.penalty = 1000;
save('LIVER.mat');

load('PROSTATE.mat');
save('PROSTATE_copy.mat');
    cst{1,6} = {y};
        cst{1,6}{1,1}.parameters = {50};
        cst{1,6}{1,1}.penalty = 300;
    cst{6,6} = {x};
        cst{6,6}{1,1}.parameters = {68};
        cst{6,6}{1,1}.penalty = 1000;
    cst{7,6} = {x};
        cst{7,6}{1,1}.parameters = {56};
        cst{7,6}{1,1}.penalty = 1000;
    cst{8,6} = {y};
        cst{8,6}{1,1}.parameters = {50};
        cst{8,6}{1,1}.penalty = 300;
    cst{9,6} = {y};
        cst{9,6}{1,1}.parameters = {30};
        cst{9,6}{1,1}.penalty = 100;
save('PROSTATE.mat');

load('TG119.mat');
save('TG119_copy.mat');
    cst{1,6} = {y};
        cst{1,6}{1,1}.parameters = {25};
        cst{1,6}{1,1}.penalty = 300;
    cst{2,6} = {x};
        cst{2,6}{1,1}.parameters = {50};
        cst{2,6}{1,1}.penalty = 1000;
    cst{3,6} = {y};
        cst{3,6}{1,1}.parameters = {30};
        cst{3,6}{1,1}.penalty = 100;
save('TG119.mat');
    
    
    
    
 