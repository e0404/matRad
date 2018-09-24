function report = matRad_unitTest_sizeCheck(verVars, cst, ct, stf,pln, dij, resultGUI)

report.status = 1;

status = [];
status = [status; nish_comparevars(cst, verVars.cst)];
status = [status; nish_comparevars(ct, verVars.ct)];
status = [status; nish_comparevars(stf, verVars.stf)];
status = [status; nish_comparevars(pln, verVars.pln)];
status = [status; nish_comparevars(dij, verVars.dij)];
status = [status; nish_comparevars(resultGUI, verVars.resultGUI)];

report.status = all(status);


end

