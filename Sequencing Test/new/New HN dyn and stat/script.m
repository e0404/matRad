dij = matRad_calcPhotonDose(ct,stf,pln,cst);

resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf,0);

pln.dynamic = true;
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
save('dynamic','resultGUI')

pln.dynamic = false;
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln,0);
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,stf);
save('static','resultGUI')

recalc.pln = pln;
recalc.dynamic = true;
recalc.interpNew = true;
recalc.dijNew = true;
dir = pwd;
cd('C:\Users\eric\Documents\GitHub\matRad\Resolution Test\H&N')
recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo,0,dij);
cd(dir)
save('static recalc','recalc')


