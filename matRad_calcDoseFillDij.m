% save computation time and memory 
% by sequentially filling the sparse matrix dose.dij from the cell array
if mod(counter,numOfBixelsContainer) == 0 || counter == dij.totalNumOfBixels

    if calcDoseDirect
        if isfield(stf(1).ray(1),'weight') && numel(stf(i).ray(j).weight) >= k

            % score physical dose
            dij.physicalDose{1}(:,i) = dij.physicalDose{1}(:,i) + stf(i).ray(j).weight(k) * doseTmpContainer{1,1};

            if isfield(dij,'mLETDose')
                dij.mLETDose{1}(:,i) = dij.mLETDose{1}(:,i) + stf(i).ray(j).weight(k) * letDoseTmpContainer{1,1};
            end

            if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
                && strcmp(pln.radiationMode,'carbon')

                % score alpha and beta matrices
                dij.mAlphaDose{1}(:,i)    = dij.mAlphaDose{1}(:,i) + stf(i).ray(j).weight(k) * alphaDoseTmpContainer{1,1};
                dij.mSqrtBetaDose{1}(:,i) = dij.mSqrtBetaDose{1}(:,i) + stf(i).ray(j).weight(k) * betaDoseTmpContainer{1,1};

            end
        else

            error(['No weight available for beam ' num2str(i) ', ray ' num2str(j) ', bixel ' num2str(k)]);

        end
    else

        dij.physicalDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [doseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];

        if isfield(dij,'mLETDose')
            dij.mLETDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [letDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
        end

        if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
            && strcmp(pln.radiationMode,'carbon')

            dij.mAlphaDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [alphaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
            dij.mSqrtBetaDose{1}(:,(ceil(counter/numOfBixelsContainer)-1)*numOfBixelsContainer+1:counter) = [betaDoseTmpContainer{1:mod(counter-1,numOfBixelsContainer)+1,1}];
        end

    end
    
end