function tk_drawShapes(shapeInfo,pln)

% xMax = max(shapeInfo.x_r)+2;
% xMin = min(shapeInfo.x_l)-2;
xMax = max(abs([shapeInfo.x_l; shapeInfo.x_r]));
xMin = 0;
numOfBeams=pln.numOfBeams;
shapeOffset=0;
offset=1;
% shapeInfo.zPos = -shapeInfo.zPos;

for m=1:pln.numOfBeams
    
    figure

    for l=1:shapeInfo.numOfShapes(m)

        numOfLines = numel(shapeInfo.shapeIx(shapeInfo.shapeIx == l));

        subplot(floor(shapeInfo.numOfShapes(m)/2),... 
            ceil(shapeInfo.numOfShapes(m)/2),l)
%         grid on
%         set(gca,'XTick',0:1:xMax)
%         axis([xMin xMax min(shapeInfo.zPos)-pln.bixelWidth/2 ...
%             max(shapeInfo.zPos)+pln.bixelWidth/2])
%         xlim([xMin xMax])
        xlim([0 xMax])
%         for i=offset:offset+numOfLines-1
%             hold on
%             fill([xMin shapeInfo.x_l(i) shapeInfo.x_l(i) xMin],...
%                 [shapeInfo.zPos(i)-pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)-pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)+pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)+pln.bixelWidth/2],'b')
%             fill([shapeInfo.x_r(i) xMax xMax shapeInfo.x_r(i)],...
%                 [shapeInfo.zPos(i)-pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)-pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)+pln.bixelWidth/2 ...
%                 shapeInfo.zPos(i)+pln.bixelWidth/2],'b')
%         end
        for i=offset:offset+numOfLines-1
            hold on
            fill([xMin shapeInfo.x_l(i) shapeInfo.x_l(i) xMin],...
                [i-1/2 ...
                i-1/2 ...
                i+1/2 ...
                i+1/2],'b')
            fill([shapeInfo.x_r(i) xMax xMax shapeInfo.x_r(i)],...
                [i-1/2 ...
                i-1/2 ...
                i+1/2 ...
                i+1/2],'b')
        end

        offset = offset + numOfLines;   
    end

    shapeOffset = shapeOffset + max(shapeInfo.numOfShapes);
end

end