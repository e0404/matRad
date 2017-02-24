function multScen = matRad_setMultScen(multScen)

% set shift scenarios
if isequal(multScen.shiftGenType,'equidistant')
    if isequal(multScen.shiftCombType,'individual')
        if isequal(multScen.shiftGen1DIsotropy,'+-')

            if sum(mod(multScen.numOfShiftScen,2)) > 0
                warning('number of shift scenarios must be a even numer due to symmetricity - affected dimension will be decreased by 1 scenario');
                multScen.numOfShiftScen = multScen.numOfShiftScen - mod(multScen.numOfShiftScen,2);
            end
            
            deltaShift = multScen.shiftSize./(multScen.numOfShiftScen./2);
            multScen.shifts  = [];

            for i = 1:3 % loop over x,y,z dimension 
                shifts      = zeros(3,multScen.numOfShiftScen(i));
                shiftstmp   = -multScen.shiftSize(i):deltaShift(i):multScen.shiftSize(i);
                shiftstmp   = shiftstmp(shiftstmp ~= 0);
                shifts(i,:) = shiftstmp;

                multScen.shifts   = [multScen.shifts, shifts];
            end
            multScen.shifts = [zeros(3,1), multScen.shifts]; % add nominal scenario

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
        if multScen.numOfShiftScen(1) == multScen.numOfShiftScen(2) &&...
           multScen.numOfShiftScen(2) == multScen.numOfShiftScen(3)
       
            if sum(mod(multScen.numOfShiftScen,2)) > 0
                warning('number of shift scenarios must be a even numer due to symmetricity - affected dimension will be decreased by 1 scenario');
                multScen.numOfShiftScen = multScen.numOfShiftScen - mod(multScen.numOfShiftScen,2);
            end
            
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
            
        else
            
            error('chose same number of shifts in x,y and z if shift comb type equals "combined"')
            
        end  
        
    elseif isequal(multScen.shiftCombType,'allcombined')

        deltaShift = multScen.shiftSize./(multScen.numOfShiftScen./2);
        multScen.shifts  = [];

        shiftsx = -multScen.shiftSize(1):deltaShift(1):multScen.shiftSize(1);
        %shiftsx = shiftsx(shiftsx ~= 0);
        shiftsy = -multScen.shiftSize(2):deltaShift(2):multScen.shiftSize(2);
        %shiftsy = shiftsy(shiftsy ~= 0);
        shiftsz = -multScen.shiftSize(3):deltaShift(3):multScen.shiftSize(3);
        %shiftsz = shiftsz(shiftsz ~= 0);

        multScen.shifts = [];
        
        for i = 1:length(shiftsx)
            for j = 1:length(shiftsy)
                for k = 1:length(shiftsz)
                    if shiftsx(i) == 0 && shiftsy(j) == 0 && shiftsz(k) == 0
                        multScen.shifts = multScen.shifts;
                    else
                        multScen.shifts = [multScen.shifts,[shiftsx(i);shiftsy(j);shiftsz(k)]];
                    end
                end
            end
        end
        
        multScen.shifts = [zeros(3,1), multScen.shifts];
        
    end
    
elseif isequal(multScen.shiftGenType,'sampled')
    
    rng(0);
    multScen.shifts  = [];
    
    if isequal(multScen.shiftCombType,'individual')
        if isequal(multScen.shiftGen1DIsotropy,'+-')

            for i = 1:3
                shifts      = zeros(3,multScen.numOfShiftScen(i));
                shifts(i,:) = multScen.shiftSD(i).*randn(1,multScen.numOfShiftScen(i));

                multScen.shifts   = [multScen.shifts, shifts];
            end
            
            multScen.shifts = [zeros(3,1), multScen.shifts];
            
        elseif isequal(multScen.shiftGen1DIsotropy,'+') ||...
               isequal(multScen.shiftGen1DIsotropy,'-')
           
           error('"+" or "-" 1D isotropy not supported in case of sampled shifts')
           
        end
    elseif isequal(multScen.shiftCombType,'combined')
        if multScen.numOfShiftScen(1) == multScen.numOfShiftScen(2) &&...
           multScen.numOfShiftScen(2) == multScen.numOfShiftScen(3)
                
            if isequal(multScen.shiftGen1DIsotropy,'+-')
                
                shiftsx = multScen.shiftSD(1).*randn(1,multScen.numOfShiftScen(1));
                shiftsy = multScen.shiftSD(2).*randn(1,multScen.numOfShiftScen(2));
                shiftsz = multScen.shiftSD(3).*randn(1,multScen.numOfShiftScen(3));
                
                multScen.shifts = [shiftsx; shiftsy; shiftsz];
                multScen.shifts = [zeros(3,1), multScen.shifts];

            elseif isequal(multScen.shiftGen1DIsotropy,'+') ||...
                   isequal(multScen.shiftGen1DIsotropy,'-')

               error('"+" or "-" 1D isotropy not supported in case of sampled shifts')

            end
        else
            
            error('chose same number of shifts in x,y and z if shift comb type equals "combined"')
            
        end
        
    end
   
end



% set range scenarios
if multScen.numOfRangeShiftScen > 0
    if multScen.absRangeShift == 0
        multScen.absRangeShifts = zeros(1,multScen.numOfRangeShiftScen+1);
    else
        deltaAbsRangeShift      = multScen.absRangeShift/(multScen.numOfRangeShiftScen/2);
        multScen.absRangeShifts = -multScen.absRangeShift:deltaAbsRangeShift:multScen.absRangeShift;
        multScen.absRangeShifts = multScen.absRangeShifts(multScen.absRangeShifts ~= 0);
        multScen.absRangeShifts = [0, multScen.absRangeShifts];
    end

    if multScen.relRangeShift == 0
        multScen.relRangeShifts = zeros(1,multScen.numOfRangeShiftScen+1);
    else
        deltaRelRangeShift      = multScen.relRangeShift/(multScen.numOfRangeShiftScen/2);
        multScen.relRangeShifts = -multScen.relRangeShift:deltaRelRangeShift:multScen.relRangeShift;
        multScen.relRangeShifts = multScen.relRangeShifts(multScen.relRangeShifts ~= 0);
        multScen.relRangeShifts = [0, multScen.relRangeShifts/100];
    end
    
elseif multScen.numOfRangeShiftScen == 0
    multScen.absRangeShifts = 0;
    multScen.relRangeShifts = 0;
end

% set total number of shift scnarios
multScen.numOfShiftScen      = size(multScen.shifts,2);
multScen.numOfRangeShiftScen = numel(multScen.absRangeShifts);
multScen.numOfScen           = multScen.numOfShiftScen + multScen.numOfRangeShiftScen -1; %substract nominal scenario of range errors

vMu = 0; vSD = 0; mShifts = 0;

% calculate probabilities of scenarios
if isequal(multScen.shiftGenType,'sampled')
    multScen.shiftScenProb = repmat(1/multScen.numOfScen,1,multScen.numOfScen);
elseif isequal(multScen.shiftGenType,'equidistant')
   
   if multScen.numOfShiftScen > 1
        vMu = [0 0 0];
        vSD = multScen.shiftSD;
        mShifts = multScen.shifts;   
    else
        mShifts = [0 0 0];
    end
    
    if multScen.numOfRangeShiftScen > 1
       rangeScen = multScen.numOfRangeShiftScen-1;
       mShifts(:,end+1:end+rangeScen) = 0;
       
       if sum(multScen.absRangeShifts ~= 0) > 1
           vMu = [vMu 0]; vSD = [vSD multScen.rangeAbsSD];
           absRangeShifts = multScen.absRangeShifts(multScen.absRangeShifts~= 0);
           for i = 1:2:numel(absRangeShifts)
              mShifts(end+1,end-rangeScen + i) = absRangeShifts(i);
              mShifts(end,end-rangeScen + i + 1)   = absRangeShifts(i+1);
           end
       end
       
       if sum(multScen.relRangeShifts ~= 0) > 1
           vMu = [vMu 0];  vSD = [vSD multScen.rangeRelSD/100];
           relRangeShifts = multScen.relRangeShifts(multScen.relRangeShifts~= 0);
           for i = 1:2:numel(relRangeShifts)
              mShifts(end+1,end-rangeScen + i) = relRangeShifts(i);
              mShifts(end,end-rangeScen + i + 1)   = relRangeShifts(i+1);
           end
       end
       
    else
       multScen.rangeScenProb = 1;
    end
    
    
    multScen.ScenProb  = matRad_calcScenProb(vMu,vSD,mShifts,'probBins','normDist');
       
else
   error(['multScen.shiftGenType :' multScen.shiftGenType ' is not implemented!']);
end



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
multScen.totalNumOfScen = sum(multScen.ScenCombMask(:));

end








