function shapeInfo = hackOpenShapes(shapeInfo)

for i = 1:size(shapeInfo.beam,2)
    
    shapeInfo.beam(i).numOfShapes = 1;
    
    shapeInfo.beam(i).shape(2:end) = [];
    
    shapeInfo.beam(i).shape.leftLeafPos  = shapeInfo.beam(i).lim_l;
    shapeInfo.beam(i).shape.rightLeafPos = shapeInfo.beam(i).lim_r;
    
    if i == 1
        shapeInfo.beam(i).shape.index = 1 + size(shapeInfo.beam,2);
    else
        shapeInfo.beam(i).shape.index = shapeInfo.beam(i-1).shape.index + numel(shapeInfo.beam(i-1).shape.leftLeafPos);
    end
    
end

shapeInfoVec = tk_shapeInfo2Vect(shapeInfo);

numOfLeafPara = (numel(shapeInfoVec)-2)/2;

shapeInfoVec(size(shapeInfo.beam,2)+[1:numOfLeafPara]) = shapeInfoVec(size(shapeInfo.beam,2)+[1:numOfLeafPara]) + shapeInfo.bixelWidth/10;
shapeInfoVec(size(shapeInfo.beam,2)+1+numOfLeafPara:end) = shapeInfoVec(size(shapeInfo.beam,2)+1+numOfLeafPara:end) - shapeInfo.bixelWidth/10;

shapeInfo = tk_updateShapeInfo(shapeInfo,shapeInfoVec);

end

