xMax=max(x_r);
numOfBeams=2;
shapeOffset=0;
m=1;
offset=1;

for m=1:pln.numOfBeams

    for l=1:numOfShapes(m)

        numOfLines = numel(shapeIx(shapeIx == l));

        subplot(numOfBeams,max(numOfShapes),shapeOffset+l)
        for i=offset:offset+numOfLines-1
            hold on
            fill([0 x_l(i) x_l(i) 0],[zPos(i)-pln.bixelWidth/2 zPos(i)-pln.bixelWidth/2 ...
                zPos(i)+pln.bixelWidth/2 zPos(i)+pln.bixelWidth/2],'r')
            fill([x_r(i) xMax xMax x_r(i)],[zPos(i)-pln.bixelWidth/2 zPos(i)-pln.bixelWidth/2 ...
                zPos(i)+pln.bixelWidth/2 zPos(i)+pln.bixelWidth/2],'r')
        end

        offset = offset + numOfLines;   
    end

    shapeOffset = shapeOffset + max(numOfShapes);
end