% Modify energy and RBE values here
energyIntRBE_Steward = ...
    [1e-9 1e-3 2e-3 5e-3 7.5e-3 1e-2 2e-2 3e-2 5e-2 1e-1 2e-1 3e-1 4e-1 5e-1 6e-1 7e-1 8e-1 9e-1 9.5e-1 1 2 2.5 3 3.5 4 5 10 20 30];

IntRBE_Steward = ...
    [3.375 3.375 3.375 3.37 3.37 3.365 3.35 3.34 3.32 3.28 3.21 3.1 3.08 2.98 2.92 2.86 2.81 2.79 2.82 2.88 2.58 2.53 2.59 2.86 2.85 2.85 2.85 2.85 2.85];

% Load energy intervals defined by tabulated KERMA values
load('energyIntTally.mat')

% Inerpolate RBE values for given energy intervals 
IntRBE_Steward_interpolatedData = interp1(energyIntRBE_Steward,IntRBE_Steward,energyIntTally);

figure 
plot(energyIntRBE_Steward, IntRBE_Steward)
hold on
plot(energyIntTally,IntRBE_Steward_interpolatedData)

data = [energyIntTally; IntRBE_Steward_interpolatedData];

fileID = fopen('neutronRBE_StewardEtAl.txt', 'w');
fprintf(fileID,'%0.5e %0.5e\n', data)
fclose(fileID)


