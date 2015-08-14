function [f, g] = matRad_daoObjFunc(shapeInfoVect,shapeInfo,dij,cst)

% update shapeInfoVect und indVect
[shapeInfo,w,indVect] = matRad_vec2ShapeInfo(shapeInfo,shapeInfoVect);

if nargout > 1
    [f, g] = matRad_objFunc(w,dij,cst);
    
    g = matRad_getGradients(g,shapeInfo.mappingMx,indVect,shapeInfo);
    
else    
    f = matRad_objFunc(w,dij,cst);
end