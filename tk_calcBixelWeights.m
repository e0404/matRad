function w = tk_calcBixelWeights(shapeInfo,dij)
% function to calculate the bixel weights from the (updated) shapeInfo
% struct using the shape weights


totalNumOfBixels = shapeInfo.totalNumOfBixels;
beamNumVect = dij.beamNum;
w = zeros(totalNumOfBixels,1);
for i = 1:totalNumOfBixels    
    % find the bixelposition corresponding to the bixel i    
        % 1. find the corresponding beam
        beamNum = beamNumVect(i);        
        % 2. find index in the MLC map
        MLCPosInd = shapeInfo.beam(beamNum).bixelIndMap == i;    
    % add up the fluence from every shape of this beam        
        for j=1:shapeInfo.beam(beamNum).numOfShapes
            w(i) = w(i) + shapeInfo.beam(beamNum).shape(j).weight * ...
                    shapeInfo.beam(beamNum).shape(j).shapeMap(MLCPosInd);
        end    
end 

end