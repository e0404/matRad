fname = input('Dij file name?','s');

cd dij

load([fname,'_dat'])

dij = dijDat;

startInd = 1;
stopInd = 0;
for i = 1:dij.numOfBeams
    load([fname,sprintf('_beam%i',i)])
    
    numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
    stopInd = stopInd+numOfRaysPerBeam;
    
    dij.physicalDose{1}(:,startInd:stopInd) = dijDos;
    
    startInd = startInd+numOfRaysPerBeam;
end

cd ..

clear fname dijDat dijDos startInd stopInd numOfRaysPerBeam i