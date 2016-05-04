function multScen = matRad_setMultScen(multScen)

% set shift scenarios
if isequal(multScen.shiftGenType,'equidistant')
    if isequal(multScen.shiftCombType,'individual')
        if isequal(multScen.shiftGen1DIsotropy,'+-')

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

        elseif isequal(multScen.shiftGen1DIsotropy,'+')

            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen);
            multScen.shifts  = [];

            for i = 1:3
                shifts      = zeros(3,multScen.numOfShiftScen(i));
                shiftstmp   = deltaShift(i):deltaShift(i):multScen.shiftSize(i);
                shifts(i,:) = shiftstmp;

                multScen.shifts   = [multScen.shifts, shifts];
            end
            multScen.shifts = [zeros(3,1), multScen.shifts];

        elseif isequal(multScen.shiftGen1DIsotropy,'-')

            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen);
            multScen.shifts  = [];

            for i = 1:3
                shifts      = zeros(3,multScen.numOfShiftScen(i));
                shiftstmp   = -multScen.shiftSize(i):deltaShift(i):-deltaShift(i);
                shifts(i,:) = shiftstmp;

                multScen.shifts   = [multScen.shifts, shifts];
            end
            multScen.shifts = [zeros(3,1), multScen.shifts];

        end
        
    elseif isequal(multScen.shiftCombType,'combined')
        if isequal(multScen.shiftGen1DIsotropy,'+-')

            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen./2);
            multScen.shifts  = [];

            shiftsx = -multScen.shiftSize(1):deltaShift(1):multScen.shiftSize(1);
            shiftsx = shiftsx(shiftsx ~= 0);
            shiftsy = -multScen.shiftSize(2):deltaShift(2):multScen.shiftSize(2);
            shiftsy = shiftsy(shiftsy ~= 0);
            shiftsz = -multScen.shiftSize(3):deltaShift(3):multScen.shiftSize(3);
            shiftsz = shiftsz(shiftsz ~= 0);
            
            multScen.shifts = [shiftsx; shiftsy; shiftsz];
            multScen.shifts = [zeros(3,1), multScen.shifts];

        elseif isequal(multScen.shiftGen1DIsotropy,'+')

            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen);
            multScen.shifts  = [];

            shiftsx = deltaShift(1):deltaShift(1):multScen.shiftSize(1);
            shiftsy = deltaShift(2):deltaShift(2):multScen.shiftSize(2);
            shiftsz = deltaShift(3):deltaShift(3):multScen.shiftSize(3);
            
            multScen.shifts = [shiftsx; shiftsy; shiftsz];
            multScen.shifts = [zeros(3,1), multScen.shifts];

        elseif isequal(multScen.shiftGen1DIsotropy,'-')

            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen);
            multScen.shifts  = [];
            
            shiftsx = -multScen.shiftSize(1):deltaShift(1):-deltaShift(1);
            shiftsy = -multScen.shiftSize(2):deltaShift(2):-deltaShift(2);
            shiftsz = -multScen.shiftSize(3):deltaShift(3):-deltaShift(3);

            multScen.shifts = [shiftsx; shiftsy; shiftsz];
            multScen.shifts = [zeros(3,1), multScen.shifts];

        end        
    end
    
    multScen.numOfShiftScen      = max(multScen.numOfShiftScen) + 1;
    
elseif isequal(multScen.shiftGenType,'sampled')
    rng(0);
    
    multScen.shifts  = [];
    
    for i = 1:3
        shifts      = zeros(3,multScen.numOfShiftScen(i));
        %shifts(i,:) = normrnd(0,multScen.shiftSize(i),1,multScen.numOfShiftScen(i));
        shifts(i,:) = multScen.shiftSize(i).*randn(1,multScen.numOfShiftScen(i));

        multScen.shifts   = [multScen.shifts, shifts];
    end
    multScen.shifts = [zeros(3,1), multScen.shifts];
    
    multScen.numOfShiftScen      = sum(multScen.numOfShiftScen)+1;
    
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

% set totalNumOfScen
multScen.totalNumOfScen = sum(multScen.ScenCombMask);

end








