fname = input('Dij file name?','s');

cd dij

load([fname,'_dat'])

dij = dijDat;

startInd = 1;
stopInd = 0;

fprintf('\nLoading dij matrices from file.\n\n');


for i = 1:dij.numOfBeams
    load([fname,sprintf('_beam%i',i)])
    
    numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
    stopInd = stopInd+numOfRaysPerBeam;
    
    dij.physicalDose{1}(:,startInd:stopInd) = dijDos;
    
    startInd = startInd+numOfRaysPerBeam;
    
    matRad_progress(i,dij.numOfBeams);
end

cd ..

clear fname dijDat dijDos startInd stopInd numOfRaysPerBeam i