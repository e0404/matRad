if strcmp(pln.radiationMode, 'brachy')
    brachyStfGen = matRad_brachyStfGenerator(pln);
    stf = brachyStfGen.generate(ct, cst, 1);
elseif strcmp(pln.radiationMode, 'photons')
    photonStfGen = matRad_photonStfGenerator(pln);
    stf = photonStfGen.generate(ct, cst, 1);
elseif any(strcmp(pln.radiationMode, {'protons', 'carbon', 'helium'}))
    ionStfGen = matRad_ionStfGenerator(pln);
    stf = ionStfGen.generate(ct, cst, 1);
end
