fname = input('Dij file name?','s');

cd dij

dijDat = rmfield(dij,'physicalDose');
save([fname,'_dat'],'dijDat');

startInd = 1;
stopInd = 0;

fprintf('\n\nSaving dij matrices to file.\n\n');

for i = 1:dij.numOfBeams
    numOfRaysPerBeam = dij.numOfRaysPerBeam(i);
    stopInd = stopInd+numOfRaysPerBeam;
    dijDos = dij.physicalDose{1}(:,startInd:stopInd);
    startInd = startInd+numOfRaysPerBeam;
    
    save([fname,sprintf('_beam%i',i)],'dijDos');
    
    matRad_progress(i,dij.numOfBeams);
end

cd ..

clear fname dijDat dijDos startInd stopInd numOfRaysPerBeam i