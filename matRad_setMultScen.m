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

% set total number of shift scnarios
multScen.numOfShiftScen = size(multScen.shifts,2);

% calculate probabilities of single shift scenarios
if isequal(multScen.shiftGenType,'sampled')
    multScen.shiftScenProb = repmat(1/multScen.numOfShiftScen,1,multScen.numOfShiftScen);
elseif isequal(multScen.shiftGenType,'equidistant')
    if multScen.numOfShiftScen > 1
        vMu = [0 0 0];
        multScen.shiftScenProb = matRad_calcScenProb(vMu,multScen.shiftSD,multScen.shifts,'probBins','normDist');
    else
        multScen.shiftScenProb = 1;
    end
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

multScen.numOfRangeShiftScen = numel(multScen.absRangeShifts);



% calculate probabilities of range shift scenarios
if isequal(multScen.shiftGenType,'sampled')
    multScen.shiftScenProb = repmat(1/multScen.numOfShiftScen,1,multScen.numOfShiftScen);
elseif isequal(multScen.shiftGenType,'equidistant')
    if multScen.numOfShiftScen > 1
        Mu  = [0];
        SD  = [multScen.rangeRelSD]/100;
        multScen.rangeScenProb = matRad_calcScenProb(Mu,SD,multScen.relRangeShifts,'probBins','normDist',[-multScen.maxRelRangeShift +multScen.maxRelRangeShift]);
    else
        multScen.rangeScenProb = 1;
    end
end


% combine propabilites 
ScenProb          = [multScen.shiftScenProb  multScen.rangeScenProb(2:end)];
multScen.ScenProb = ScenProb/sum(ScenProb);

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








