function [cst,pln] = matRad_setPlanUncertainties(ct,cst,pln)


pln.robOpt = false;

for i = 1:size(cst,1)
  for j = 1:numel(cst{i,6})
      if ~strcmp(cst{i,6}(j).robustness,'none')
         pln.robOpt = true;
      end
  end
end
       
 %% create multiple scenario struc
if pln.robOpt
   multScen.numOfCtScen          = ct.numOfCtScen; % number of imported ct scenarios

   multScen.numOfIntSegShiftScen = 0; %1000;       % number of internal segmentation shift scnearios     

   multScen.numOfShiftScen       = [0 0 0];        % number of shifts in x y and z direction       
   multScen.shiftSize            = [3 3 3];        % maximum shift [mm]  % prostate cases 5mm otherwise 3mm
   multScen.shiftSD              = [3 3 3];        % SD of normal distribution [mm]
   multScen.shiftGenType         = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.shiftCombType        = 'individual';   % individual: no combination of shift scenarios, combined: combine shift scenarios, allcombined: create every possible shift combination
   multScen.shiftGen1DIsotropy   = '+-';           % for equidistant shifts: '+-': positive and negative, '-': negative, '+': positive shift generation 

   multScen.numOfRangeShiftScen  = 0;              % number of absolute and/or relative range scnearios. if absolute and relative range scenarios are defined then relative and absolute errors are combined
   multScen.absRangeShift        = 0;              % maximum absolute over and undershoot in mm
   multScen.rangeAbsSD           = 1;              % SD of normal distribution
   multScen.relRangeShift        = 3.5;            % maximum relative over and undershoot in %
   multScen.rangeRelSD           = 3.5;            % SD of normal distribution
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
   multScen.rangeGenType         = 'equidistant';  % equidistant: equidistant shifts, sampled: sample shifts from normal distribution
   multScen.ScenCombType         = 'individual';   % individual: no combination of scenarios, allcombined: combine all scenarios
 
end


pln.multScen = matRad_setMultScen(multScen);



%% coverage based cst manipulation
cst = matRad_robOptCstManipulation(ct,cst,pln,multScen,0,0);





end