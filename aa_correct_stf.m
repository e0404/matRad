stf2 = stf;
stf2.ray(1) = stf2.ray(124);
stf2.ray(2:end) = [];
stf2.numOfBixelsPerRay(1) = stf2.numOfBixelsPerRay(124);
stf2.numOfBixelsPerRay(2:end) = [];
stf2.numOfRays = 1;
stf2.totalNumOfBixels = 1;
stf2.numOfBixelsPerRay = 1;
stf2.ray.energy(1) = stf2.ray.energy(end);
stf2.ray.energy(2:end) = [];
stf2.ray.focusIx(1) = stf2.ray.focusIx(end);
stf2.ray.focusIx(2:end) = [];
stf2.numOfRays = 1;

resultGUI2 = resultGUI;
resultGUI2.w = 1;

resultGUI3 = resultGUI;
resultGUI3.w(56:end) = []; 