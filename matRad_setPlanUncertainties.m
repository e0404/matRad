function [cst,pln] = matRad_setPlanUncertainties(ct,cst,pln)

 %% create multiple scenario struc
if ~strcmp(pln.robOpt,'none')
   multScen.numOfCtScen          = ct.numOfCtScen; % number of imported ct scenarios

   multScen.numOfIntSegShiftScen = 0; %1000;       % number of internal segmentation shift scnearios     

   multScen.numOfShiftScen       = [0 0 0];        % number of shifts in x y and z direction       
   multScen.shiftSize            = [4 4 4];        % maximum shift [mm]
   multScen.shiftSD              = [3 3 3];        % SD of normal distribution [mm]
   multScen.shiftGenType         = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.shiftCombType        = 'individual';   % individual: no combination of shift scenarios, combined: combine shift scenarios, allcombined: create every possible shift combination
   multScen.shiftGen1DIsotropy   = '+-';           % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

   multScen.numOfRangeShiftScen  = 2;              % number of absolute and/or relative range scnearios. if absolute and relative range scenarios are defined then relative and absolute errors are combined
   multScen.maxAbsRangeShift     = 0;              % maximum absolute over and undershoot in mm
   multScen.maxRelRangeShift     = 3.5;              % maximum relative over and undershoot in %

   multScen.ScenCombType         = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
   
else
    
   multScen.numOfCtScen          = ct.numOfCtScen; % number of imported ct scenarios

   multScen.numOfIntSegShiftScen = 0; %1000;       % number of internal segmentation shift scnearios     

   multScen.numOfShiftScen       = [0 0 0];        % number of shifts in x y and z direction       
   multScen.shiftSize            = [0 0 0];        % maximum shift [mm]
   multScen.shiftSD              = [0 0 0];        % SD of normal distribution [mm]
   multScen.shiftGenType         = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.shiftCombType        = 'individual';   % individual: no combination of shift scenarios, combined: combine shift scenarios, allcombined: create every possible shift combination
   multScen.shiftGen1DIsotropy   = '+-';           % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

   multScen.numOfRangeShiftScen  = 0;              % number of absolute and/or relative range scnearios. if absolute and relative range scenarios are defined then relative and absolute errors are combined
   multScen.maxAbsRangeShift     = 0;              % maximum absolute over and undershoot in mm
   multScen.maxRelRangeShift     = 0;              % maximum relative over and undershoot in %

   multScen.ScenCombType         = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
 
end


pln.multScen                      = matRad_setMultScen(multScen);



%% coverage based cst manipulation
cst = matRad_robOptCstManipulation(ct,cst,pln,multScen,0,0);





end