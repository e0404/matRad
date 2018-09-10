function report = matRad_unitTest_sizeCheck(cst, ct, stf, dij, resultGui)

report.status = 1;
report.cst.size2 = size(cst, 2);

if(report.cst.size2 != 5)
    report.status = 0;
    warning ("Wrong size: expected 5, but got %d",report.cst.size2)
endif


end