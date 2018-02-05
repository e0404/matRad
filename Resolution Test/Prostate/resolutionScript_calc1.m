tic

angularResS = [2 4];

for angularRes = angularResS
    %for each angular resolution, proceed from the best approximation to
    %the worst
    recalc.pln = pln;
    recalc.pln.minGantryAngleRes = angularRes;
    
    
    %first time, do interpolation and dynamic fluence calculation
    fname = sprintf('%.1f degrees, dyn + interp.mat',angularRes);
    recalc.dynamic = true;
    recalc.interpNew = true;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    
    
    
    %next, do dynamic fluence and interpolation, but using old dij matrices
    fname = sprintf('%.1f degrees, dyn + interp oldDij.mat',angularRes);
    recalc.dynamic = true;
    recalc.interpNew = true;
    recalc.dijNew = false;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    
    
    %{
    %NOT SURE IT MAKES SENSE TO EVER DO THIS
    %next, do dynamic fluence but no interpolation
    fname = sprintf('%.1f degrees, dyn + Ninterp.mat',angularRes);
    recalc.dynamic = true;
    recalc.interpNew = false;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    %}
    
    
    %next, do interpolation but no dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + interp.mat',angularRes);
    recalc.dynamic = false;
    recalc.interpNew = true;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    
    
    %finally, do neither interpolation nor dynamic fluence
    fname = sprintf('%.1f degrees, Ndyn + Ninterp.mat',angularRes);
    recalc.dynamic = false;
    recalc.interpNew = false;
    recalc.dijNew = true;
    
    recalc = matRad_doseRecalc(cst,pln,recalc,ct,resultGUI.apertureInfo);
    save(fname,'resultGUI','recalc');
    
    
end

toc