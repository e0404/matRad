tic
resultGUI = matRad_fluenceOptimization(dij,cst,pln,stf,0);
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,pln.numLevels,0,pln.VMAT,0,pln);
resultGUI = matRad_directApertureOptimization(dij,cst,resultGUI.apertureInfo,resultGUI,pln,1,1);
toc