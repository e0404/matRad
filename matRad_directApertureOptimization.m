function optResult = matRad_directApertureOptimization(dij,cst,shapeInfo,visBool,varargin)

% generate initial vector for optimization
initShapeInfoVec = shapeInfo.vector;

% set objective function
objFunc =  @(x) matRad_daoObjFunc(x,shapeInfo,dij,cst);

% set projection function
projFunc = @(x) matRad_daoProjectionFunction(x,shapeInfo.limMx,shapeInfo.totalNumOfShapes);

% verify gradients
%matRad_verifyGradient(objFunc,initShapeInfoVec);

% minimize objetive function
optShapeInfoVec = matRad_projectedLBFGS(objFunc,projFunc,initShapeInfoVec,visBool,varargin);


%% calc bixel weights
% update the shapeInfoStruct
[optResult.shapeInfo,optResult.w] = matRad_vec2ShapeInfo(shapeInfo,optShapeInfoVec);

% calc dose and reshape from 1D vector to 2D array
optResult.physicalDose = reshape(dij.physicalDose*optResult.w,dij.dimensions);



