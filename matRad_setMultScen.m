function ms = matRad_setMultScen(ms)

% set shift scenarios
deltaShift = ms.maxShift./(ms.numOfShiftScen./2);
ms.shifts  = [];

for i = 1:3
    shifts      = zeros(3,ms.numOfShiftScen(i));
    shiftstmp   = -ms.maxShift(i):deltaShift(i):ms.maxShift(i);
    shiftstmp   = shiftstmp(shiftstmp ~= 0);
    shifts(i,:) = shiftstmp;
    
    ms.shifts   = [ms.shifts, shifts];
end
ms.shifts = [zeros(3,1), ms.shifts];

% set range scenarios
deltaRange = ms.maxRange/(ms.numOfRangeScen/2);
ms.range   = -ms.maxRange:deltaRange:ms.maxRange;
ms.range   = ms.range(ms.range ~= 0);
ms.range   = [0, ms.range];

% set scenario combination mask
if isequal(ms.ScenCombType,'individual')
    
    ms.ScenCombMask = false(ms.numOfCtScen, sum(ms.numOfShiftScen)+1, ms.numOfRangeScen+1);
    
    ms.ScenCombMask(:,1,1) = true; % individual ct scenarios
    ms.ScenCombMask(1,:,1) = true; % individual shift scenarios on ct scenario 1
    ms.ScenCombMask(1,1,:) = true; % individual range scenarios on ct scenario 1
    
elseif isequal(ms.ScenCombType,'allcombined')
    
    ms.ScenCombMask = true(ms.numOfCtScen, sum(ms.numOfShiftScen)+1, ms.numOfRangeScen+1);
    
end

end








