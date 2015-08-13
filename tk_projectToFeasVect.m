function [feasibleVect, isConstrActive] = tk_projectToFeasVect(shapeInfoVect,limVect,createLogBool)
% function to correct the
% shapeInfoVect[shapeWeights,LeftLeafPos,RightLeafPos] for problems:
% 1. negative weights
% 2. Overlap of left and right leaf
% 3. Leaf out of bounds

%% create logFile for the corrections
if createLogBool 
    log = fopen('logFile.txt','a');
    fprintf(log,'\n\n----------\n newEntry \n\n');
end

%%
% initialize Variables
totalNumOfShapes = sum(limVect(:,2) == inf);
totalNumOfLeafPairs = (1/2) * (size(limVect,1)-totalNumOfShapes);
isConstrActive = false(size(shapeInfoVect,1),1);
%% correct vector
% first: put all positions into limits
for i=1:length(shapeInfoVect)
    % correct negative shape weights and leaf bounds
    if shapeInfoVect(i)<=limVect(i,1)
        shapeInfoVect(i) = limVect(i,1);
        isConstrActive(i) = 1;
        if createLogBool            
            if i>totalNumOfShapes
                fprintf(log,['left_lim corrected at index ' num2str(i) '\n']);
            else
                fprintf(log,['negative weight corrected at index ' num2str(i) '\n']);
            end
        end
    elseif shapeInfoVect(i)>=limVect(i,2)
        shapeInfoVect(i) = limVect(i,2);        
        isConstrActive(i) = 1;
        if createLogBool 
            fprintf(log,['right_lim corrected at index ' num2str(i) '\n']);
        end
    end
    
end

% second: correct overlap
for i=totalNumOfShapes+1:length(shapeInfoVect)-totalNumOfLeafPairs
    if shapeInfoVect(i) >= shapeInfoVect(i+totalNumOfLeafPairs)
            overlap = shapeInfoVect(i)-shapeInfoVect(i+totalNumOfLeafPairs);
            % set leaves to the center of their overlap
            shapeInfoVect(i) = shapeInfoVect(i) - overlap/2;
            shapeInfoVect(i+totalNumOfLeafPairs) = shapeInfoVect(i);            
            isConstrActive(i) = 1;
            if createLogBool 
                fprintf(log,['overlap corrected at index ' num2str(i) ' - ' ...
                        num2str(i+totalNumOfLeafPairs) '\n']);
            end
    end
end
    

feasibleVect = shapeInfoVect;

%% close logFile
if createLogBool
    fclose(log);
end


end

