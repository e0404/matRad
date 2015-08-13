function [feasibleVect, isConstrActive] = matRad_projectToFeasVect(weightVect)
% function to correct the
% shapeInfoVect[shapeWeights,LeftLeafPos,RightLeafPos] for problems:
% 1. negative weights
%%
%% correct vector
weightVect(weightVect<0) = 0;

isConstrActive = weightVect==0;

feasibleVect = weightVect;

end

