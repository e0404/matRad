function optResult = matRad_DAOoptimization(dij,stf,cst,pln,Sequencing,visBool,varargin)

shapeInfo = tk_getParameters(Sequencing,stf,pln,0);
[shapeInfoVect, addInfoVect] = tk_shapeInfo2Vect(shapeInfo);
%[shapeInfo] = tk_updateShapeInfo(shapeInfo,shapeInfoVect);

% set objective function
objFunc =  @(x) matRad_DAOobjFunc(x,shapeInfo,addInfoVect,dij,cst);

% set projection function
limVect = tk_createLimVect(shapeInfo,addInfoVect);

% projFunc = @(x) tk_projectToFeasVect(x,limVect,1);


%% verify gradients
% matRad_verifyDAOGradient(objFunc,shapeInfoVect);

%% optimization
% create limit vector containing feasible ranges for all elements of
% shapeInfoVect
limVect = tk_createLimVect(shapeInfo,addInfoVect);
% minimize objetive function
optResult = matRad_DAOprojectedLBFGS(objFunc,shapeInfoVect,limVect,visBool,varargin);


%% calc bixel weights
% update the shapeInfoStruct
optResult.shapeInfo = tk_updateShapeInfo(shapeInfo,optResult.w);
% calculate new weight vector
optResult.w = tk_calcBixelWeights(optResult.shapeInfo,dij);

% calc dose and reshape from 1D vector to 2D array
optResult.physicalDose = reshape(dij.physicalDose*optResult.w,dij.dimensions);



