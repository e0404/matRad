% save computation time and memory 
% by sequentially filling the sparse matrix dose.dij from the cell array
if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels
    %First check if we have a weight if calcDoseDirect to not do it in the
    %innermost loop every single time
    if calcDoseDirect && (~isfield(stf(1).ray(1),'weight') || numel(stf(i).ray(j).weight) < k)
        matRad_cfg.dispError(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)]);
    end
        
    for ctScen = 1:pln.multScen.numOfCtScen
        for rangeShiftScen = 1:pln.multScen.totNumRangeScen
            if pln.multScen.scenMask(ctScen,shiftScen,rangeShiftScen)
                if calcDoseDirect
                    rayWeight = stf(i).ray(j).weight(k);
                    % score physical dose
                    dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,i) + rayWeight * doseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                    
                    if isfield(dij,'mLETDose')
                        % score LETxDose matrices
                        dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,i) + rayWeight * letDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                    end
                    
                    if pln.bioParam.bioOpt
                        % score alphaxDose and sqrt(beta)xDose matrices
                        dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,i)    = dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,i)    + rayWeight * alphaDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                        dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,i) = dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,i) + rayWeight * betaDoseTmpContainer{1,ctScen,shiftScen,rangeShiftScen};
                    end
                else
                    dijIx = (ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter;
                    containerIx = 1:mod(counter-1,numOfBixelsContainer)+1;
                    
                    % fill entire dose influence matrix
                    dij.physicalDose{ctScen,shiftScen,rangeShiftScen}(:,dijIx) = [doseTmpContainer{containerIx,ctScen,shiftScen,rangeShiftScen}];
                    
                    if isfield(dij,'mLETDose')
                        % fill entire LETxDose influence matrix
                        dij.mLETDose{ctScen,shiftScen,rangeShiftScen}(:,dijIx) = [letDoseTmpContainer{containerIx,ctScen,shiftScen,rangeShiftScen}];
                    end
                    
                    if pln.bioParam.bioOpt
                        % fill entire alphaxDose influence and sqrt(beta)xDose influence matrices
                        dij.mAlphaDose{ctScen,shiftScen,rangeShiftScen}(:,dijIx)    = [alphaDoseTmpContainer{containerIx,ctScen,shiftScen,rangeShiftScen}];
                        dij.mSqrtBetaDose{ctScen,shiftScen,rangeShiftScen}(:,dijIx) = [betaDoseTmpContainer{containerIx,ctScen,shiftScen,rangeShiftScen}];
                    end
                end
            end
        end
    end    
end