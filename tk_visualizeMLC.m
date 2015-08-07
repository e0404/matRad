function tk_visualizeMLC(shapeInfo,pln)

% define what to draw
mode = 'physical'; % options: 'physical','leafNum'

% global parameters
numOfBeams = pln.numOfBeams;
bixelWidth = pln.bixelWidth;

if strcmp(mode,'physical')
%   loop over all beams
    for i=1:numOfBeams
%       open new figure for every beam
        figure
        
%       get the MLC dimensions for this beam
        minX = shapeInfo.beam(i).MLCWindow(1);
        maxX = shapeInfo.beam(i).MLCWindow(2);        
        
%       loop over all shapes of the beam 
        for j=1:shapeInfo.beam(i).numOfShapes
            
%           creating subplots
            subplotColumns = ceil(shapeInfo.beam(i).numOfShapes/2);
            remainingPlots = shapeInfo.beam(i).numOfShapes - subplotColumns;
            if remainingPlots == 0
                subplotLines = 1;
            else
                subplotLines = 2;
            end
            subplot(subplotLines,subplotColumns,j)

            title(['Beam: ' num2str(i) ' Shape: ' num2str(j)],...
                        'Fontsize',10)
            hold on
%           loop over all active leaf pairs
            for k=1:shapeInfo.beam(i).numOfActiveLeafPairs
                fill([minX shapeInfo.beam(i).shape(j).leftLeafPos(k) ...
                    shapeInfo.beam(i).shape(j).leftLeafPos(k) minX],...
                    [shapeInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)+ bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)+ bixelWidth/2],'b')
                fill([shapeInfo.beam(i).shape(j).rightLeafPos(k) ...
                    maxX maxX ...
                    shapeInfo.beam(i).shape(j).rightLeafPos(k)],...
                    [shapeInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)- bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)+ bixelWidth/2 ...
                    shapeInfo.beam(i).leafPairPos(k)+ bixelWidth/2],'b')                
            end
            axis tight
            xlabel('horiz. pos. [mm]')
            ylabel('vert. pos. [mm]')
        end
        
    end
    
end

if strcmp(mode,'leafNum')
%   loop over all beams
    for i=1:numOfBeams
%       open new figure for every beam
        figure
        
%       get the MLC dimensions for this beam
        minX = shapeInfo.beam(i).MLCWindow(1);
        maxX = shapeInfo.beam(i).MLCWindow(2);     
        
%       get the active leaf Pairs
        activeLeafInd = find(shapeInfo.beam(i).isActiveLeafPair);
%       the leaf indices have to be flipped in order to fit to the order of
%       the leaf positions (1st row of leafPos is lowest row in physical
%       coordinates
        activeLeafInd = flipud(activeLeafInd); % flip to fit to the 

%       creating subplots
        subplotColumns = ceil(shapeInfo.beam(i).numOfShapes/2);
        remainingPlots = shapeInfo.beam(i).numOfShapes - subplotColumns;
        if remainingPlots == 0
            subplotLines = 1;
        else
            subplotLines = 2;
        end        
        
%       loop over all shapes of the beam 
        for j=1:shapeInfo.beam(i).numOfShapes
            
            subplot(subplotLines,subplotColumns,j)

            title(['Beam: ' num2str(i) ' Shape: ' num2str(j)],...
                        'Fontsize',10)
            hold on
%           loop over all active leaf pairs
            for k=1:shapeInfo.beam(i).numOfActiveLeafPairs
                fill([minX shapeInfo.beam(i).shape(j).leftLeafPos(k) ...
                    shapeInfo.beam(i).shape(j).leftLeafPos(k) minX],...
                    [activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) + 1/2 ...
                    activeLeafInd(k) + 1/2],'b') 
                fill([shapeInfo.beam(i).shape(j).rightLeafPos(k) ...
                    maxX maxX ...
                    shapeInfo.beam(i).shape(j).rightLeafPos(k)],...
                    [activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) - 1/2 ...
                    activeLeafInd(k) + 1/2 ...
                    activeLeafInd(k) + 1/2],'b')                
            end
            axis tight
            axis ij
            xlabel('horiz. pos. [mm]')
            ylabel('leaf pair # [mm]')
        end
        
    end
    
end





end