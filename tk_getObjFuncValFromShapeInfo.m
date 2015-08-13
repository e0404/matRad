function [w, objFuncVal, gradient] = tk_getObjFuncValFromShapeInfo(shapeInfo,dij,cst)

objFunc =  @(x) matRad_objFunc(x,dij,cst);

% initializing variables
totalNumOfBixels = shapeInfo.totalNumOfBixels;
beamNumVect = dij.beamNum;
w = zeros(totalNumOfBixels,1);

%% 1. Get the weight vector from the shapes of each beam
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

%% 2. Get objective function value and gradient

[objFuncVal, gradient] = objFunc(w);

%matRad_verifyGradient(objFunc,length(w));

end