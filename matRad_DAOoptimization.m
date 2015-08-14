function optResult = matRad_DAOoptimization(dij,stf,cst,pln,Sequencing,visBool,varargin)

[shapeInfo, shapeInfoVect] = matRad_getSeqParameters(Sequencing,stf,pln,0);
% [shapeInfoVect, addInfoVect] = matRad_shapeInfo2Vect(shapeInfo);
%[shapeInfo] = matRad_updateShapeInfo(shapeInfo,shapeInfoVect);

% set objective function
objFunc =  @(x) matRad_DAOobjFunc(x,shapeInfo,dij,cst);

% set projection function
% limVect = matRad_createLeafLimVect(shapeInfo,addInfoVect);

projFunc = @(x) matRad_DAOprojectToFeasVect(x,shapeInfo.limVect,1);


%% verify gradients
% matRad_verifyDAOGradient(objFunc,shapeInfoVect);

%% optimization
% create limit vector containing feasible ranges for all elements of

% minimize objetive function
optResult = matRad_DAOprojectedLBFGS(objFunc,projFunc,shapeInfoVect,visBool,varargin);


%% calc bixel weights
% update the shapeInfoStruct
optResult.shapeInfo = matRad_updateShapeInfo(shapeInfo,optResult.w);
% calculate new weight vector
optResult.w = matRad_calcBixelWeights(optResult.shapeInfo,dij);

% calc dose and reshape from 1D vector to 2D array
optResult.physicalDose = reshape(dij.physicalDose*optResult.w,dij.dimensions);



