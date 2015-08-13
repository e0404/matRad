function tk_drawShapes(shapeInfo,pln)

% xMax = max(shapeInfo.x_r)+2;
% xMin = min(shapeInfo.x_l)-2;
xMax = max(abs([shapeInfo.x_l; shapeInfo.x_r]));
xMin = 0;

% define wether bixel from 0.5 to 1.5 or from 0 to 1
% 0.5: from 0.5 to 1.5 (like imagesc)
% 0: from 0 to 1
bixelOffset = 0; 

% define what to draw (indices or physical coordinates)
drawMode = 'leafNum'; %options: 'index','physical' 'leafNum'

numOfBeams=pln.numOfBeams;
shapeOffset=0;
offset=1;
% shapeInfo.zPos = -shapeInfo.zPos;

for m=1:numOfBeams
    
    figure

    for l=1:shapeInfo.numOfShapes(m)

        numOfLines = numel(shapeInfo.shapeIx(shapeInfo.shapeIx == l+shapeOffset));
        
        % creating subplots
        % lines
        subplotLines = floor(shapeInfo.numOfShapes(m)/2);
        subplotRows = ceil(shapeInfo.numOfShapes(m)/2);
        if shapeInfo.numOfShapes(m)==3
            subplot(2,2,l)
        else
            subplot(subplotLines,... 
            subplotRows,l)
        end

%         set properties
%         grid on
%         set(gca,'XTick',0:1:xMax)
%         axis([xMin xMax min(shapeInfo.zPos)-pln.bixelWidth/2 ...
%             max(shapeInfo.zPos)+pln.bixelWidth/2])
%         xlim([xMin xMax])

%         xlim([0.5 xMax+0.5])
%         ylim([offset-0.5 offset+numOfLines-1+0.5])
        
        %physical coordinates
        if strcmp(drawMode, 'physical')
            for i=offset:offset+numOfLines-1
                hold on
                fill([xMin shapeInfo.x_l(i) shapeInfo.x_l(i) xMin],...
                    [shapeInfo.leafIx(i)-pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)-pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)+pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)+pln.bixelWidth/2],'b')
                fill([shapeInfo.x_r(i) xMax xMax shapeInfo.x_r(i)],...
                    [shapeInfo.leafIx(i)-pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)-pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)+pln.bixelWidth/2 ...
                    shapeInfo.leafIx(i)+pln.bixelWidth/2],'b')
            end
            axis ij
            axis tight
        end
        
        % using indices for visualization
        if strcmp(drawMode, 'index')

%             xMin = shapeInfo.x_min(offset);
            xMax = shapeInfo.x_max(offset);

            hold on
            for i=offset:offset+numOfLines-1
                fill([bixelOffset shapeInfo.x_l(i)+bixelOffset shapeInfo.x_l(i)+bixelOffset bixelOffset],...
                    [i-1/2 ...
                    i-1/2 ...
                    i+1/2 ...
                    i+1/2],'b')
                fill([shapeInfo.x_r(i)+bixelOffset xMax+bixelOffset xMax+bixelOffset shapeInfo.x_r(i)+bixelOffset],...
                    [i-1/2 ...
                    i-1/2 ...
                    i+1/2 ...
                    i+1/2],'b')
            end
            axis ij
            axis tight
            offset = offset + numOfLines;
        end
        
         % using leafNumbers for visualization
        if strcmp(drawMode, 'leafNum')

%             xMin = shapeInfo.x_min(offset);
            xMax = shapeInfo.x_max(offset);

            hold on
            for i=offset:offset+numOfLines-1
                fill([bixelOffset shapeInfo.x_l(i)+bixelOffset shapeInfo.x_l(i)+bixelOffset bixelOffset],...
                    [shapeInfo.leafIx(i)-1/2 ...
                    shapeInfo.leafIx(i)-1/2 ...
                    shapeInfo.leafIx(i)+1/2 ...
                    shapeInfo.leafIx(i)+1/2],'b')
                fill([shapeInfo.x_r(i)+bixelOffset xMax+bixelOffset xMax+bixelOffset shapeInfo.x_r(i)+bixelOffset],...
                    [shapeInfo.leafIx(i)-1/2 ...
                    shapeInfo.leafIx(i)-1/2 ...
                    shapeInfo.leafIx(i)+1/2 ...
                    shapeInfo.leafIx(i)+1/2],'b')
            end
            axis ij
            axis tight
            offset = offset + numOfLines;
        end
    end

    shapeOffset = shapeOffset + shapeInfo.numOfShapes(m);
end

end