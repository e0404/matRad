function multScen = matRad_setMultScen(multScen)

% set shift scenarios
if isequal(multScen.shiftGenType,'equidistant')
    deltaShift = multScen.shiftSize./(multScen.numOfShiftScen./2);
    multScen.shifts  = [];

    for i = 1:3
        shifts      = zeros(3,multScen.numOfShiftScen(i));
        shiftstmp   = -multScen.shiftSize(i):deltaShift(i):multScen.shiftSize(i);
        shiftstmp   = shiftstmp(shiftstmp ~= 0);
        shifts(i,:) = shiftstmp;

        multScen.shifts   = [multScen.shifts, shifts];
    end
    multScen.shifts = [zeros(3,1), multScen.shifts];
    
elseif isequal(multScen.shiftGenType,'sampled')
    multScen.shifts  = [];
    
    for i = 1:3
        shifts      = zeros(3,multScen.numOfShiftScen(i));
        shifts(i,:) = normrnd(0,multScen.shiftSize(i),1,multScen.numOfShiftScen(i));

        multScen.shifts   = [multScen.shifts, shifts];
    end
    multScen.shifts = [zeros(3,1), multScen.shifts];
    
end

% set range scenarios
if multScen.numOfRangeShiftScen > 0
    if multScen.maxAbsRangeShift == 0
        multScen.absRangeShifts = zeros(1,multScen.numOfRangeShiftScen+1);
    else
        deltaAbsRangeShift      = multScen.maxAbsRangeShift/(multScen.numOfRangeShiftScen/2);
        multScen.absRangeShifts = -multScen.maxAbsRangeShift:deltaAbsRangeShift:multScen.maxAbsRangeShift;
        multScen.absRangeShifts = multScen.absRangeShifts(multScen.absRangeShifts ~= 0);
        multScen.absRangeShifts = [0, multScen.absRangeShifts];
    end

    if multScen.maxRelRangeShift == 0
        multScen.relRangeShifts = zeros(1,multScen.numOfRangeShiftScen+1);
    else
        deltaRelRangeShift      = multScen.maxRelRangeShift/(multScen.numOfRangeShiftScen/2);
        multScen.relRangeShifts = -multScen.maxRelRangeShift:deltaRelRangeShift:multScen.maxRelRangeShift;
        multScen.relRangeShifts = multScen.relRangeShifts(multScen.relRangeShifts ~= 0);
        multScen.relRangeShifts = [0, multScen.relRangeShifts/100];
    end
    
elseif multScen.numOfRangeShiftScen == 0
    multScen.absRangeShifts = 0;
    multScen.relRangeShifts = 0;
end

% set correct number of scenarios
multScen.numOfShiftScen      = sum(multScen.numOfShiftScen)+1;
multScen.numOfRangeShiftScen = multScen.numOfRangeShiftScen+1;

% set scenario combination mask
if isequal(multScen.ScenCombType,'individual')
    
    multScen.ScenCombMask = false(multScen.numOfCtScen, multScen.numOfShiftScen, multScen.numOfRangeShiftScen);
    
    multScen.ScenCombMask(:,1,1) = true; % individual ct scenarios
    multScen.ScenCombMask(1,:,1) = true; % individual shift scenarios on ct scenario 1
    multScen.ScenCombMask(1,1,:) = true; % individual range scenarios on ct scenario 1
    
elseif isequal(multScen.ScenCombType,'allcombined')
    
    multScen.ScenCombMask = true(multScen.numOfCtScen, multScen.numOfShiftScen, multScen.numOfRangeShiftScen);
    
end

end








