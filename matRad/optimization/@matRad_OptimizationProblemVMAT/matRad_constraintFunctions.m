function c = matRad_constraintFunctions(optiProb,apertureInfoVec,dij,cst)
% matRad IPOPT callback: constraint function for VMAT
%
% call
%   c = matRad_daoObjFunc(apertueInfoVec,dij,cst)
%
% input
%   apertueInfoVec: aperture info vector
%   dij:            dose influence matrix
%   cst:            matRad cst struct
%   options:        option struct defining the type of optimization
%
% output
%   c:              value of constraints
%
% Reference
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DAO constraint function calculation
c_dos_dao = matRad_constraintFunctions@matRad_OptimizationProblemDAO(optiProb,apertureInfoVec,dij,cst);

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
%This update should have taken place if necessary in the call above
%if ~isequal(apertureInfoVec,optiProb.apertureInfo.apertureVector)
%    optiProb.apertureInfo = matRad_daoVec2ApertureInfo(optiProb.apertureInfo,apertureInfoVec);
%end
apertureInfo = optiProb.apertureInfo;

% We need the leaf pos
leftLeafPos  = apertureInfoVec((1:apertureInfo.totalNumOfLeafPairs)+apertureInfo.totalNumOfShapes);
rightLeafPos = apertureInfoVec((1+apertureInfo.totalNumOfLeafPairs+apertureInfo.totalNumOfShapes):(apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2));

% values of times spent in an arc surrounding the optimized angles (full
% arc/dose influence arc)
timeDAOBorderAngles = apertureInfoVec(((apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2)+1):end);
timeDoseBorderAngles = timeDAOBorderAngles.*[apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFacCurr]';

if apertureInfo.propVMAT.continuousAperture
    % Using the dynamic fluence calculation, we have the leaf positions in
    % the vector be the leaf positions at the borders of the Dij arcs (for optimized angles only).
    % Therefore we must also use the times between the borders of the Dij
    % arc (for optimized angles only).
    timeFac = [apertureInfo.propVMAT.beam.timeFac]';
    deleteInd = timeFac == 0;
    timeFac(deleteInd) = [];
    
    i = [apertureInfo.propVMAT.beam.timeFacInd]';
    i(deleteInd) = [];
    
    j = repelem(1:apertureInfo.totalNumOfShapes,1,3);
    j(deleteInd) = [];
    
    timeFacMatrix = sparse(i,j,timeFac,max(i),apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix*timeDAOBorderAngles;
    
    % prep
    leftLeafSpeed   = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
    rightLeafSpeed  = zeros(apertureInfo.propVMAT.numLeafSpeedConstraint*apertureInfo.beam(1).numOfActiveLeafPairs,1);
    
    offset      = 0;
    shapeInd    = 1;
    
    for i = 1:numel(apertureInfo.beam)
        % loop over beams
        n = apertureInfo.beam(i).numOfActiveLeafPairs;
        
        if ~isempty(apertureInfo.propVMAT.beam(i).leafConstMask)
            
            % get vector indices
            if apertureInfo.propVMAT.beam(i).DAOBeam
                % if it's a DAO beam, use own vector offset
                vectorIx_LI = apertureInfo.beam(i).shape(1).vectorOffset(1) + ((1:n)-1);
                vectorIx_LF = apertureInfo.beam(i).shape(1).vectorOffset(2) + ((1:n)-1);
            else
                % otherwise, use vector offset of previous and next
                % beams
                vectorIx_LI = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).vectorOffset(2) + ((1:n)-1);
                vectorIx_LF = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).vectorOffset(1) + ((1:n)-1);
            end
            vectorIx_RI = vectorIx_LI+apertureInfo.totalNumOfLeafPairs;
            vectorIx_RF = vectorIx_LF+apertureInfo.totalNumOfLeafPairs;
            
            % extract leaf positions, time
            leftLeafPos_I   = apertureInfoVec(vectorIx_LI);
            rightLeafPos_I  = apertureInfoVec(vectorIx_RI);
            leftLeafPos_F   = apertureInfoVec(vectorIx_LF);
            rightLeafPos_F  = apertureInfoVec(vectorIx_RF);
            t               = timeBNOptAngles(shapeInd);
            
            % determine indices
            indInConVec = offset+(1:n);
            
            % calc speeds
            leftLeafSpeed(indInConVec)      = abs(leftLeafPos_F-leftLeafPos_I)./t;
            rightLeafSpeed(indInConVec)     = abs(rightLeafPos_F-rightLeafPos_I)./t;
            
            % update offset
            offset = offset+n;
            
            % increment shapeInd only for beams which have transtion
            % defined
            shapeInd = shapeInd+1;
        end
    end
    
    c_lfspd = [leftLeafSpeed; rightLeafSpeed];
else
    
    i = sort(repmat(1:(apertureInfo.totalNumOfShapes-1),1,2));
    j = sort(repmat(1:apertureInfo.totalNumOfShapes,1,2));
    j(1) = [];
    j(end) = [];
    
    timeFac = [apertureInfo.propVMAT.beam([apertureInfo.propVMAT.beam.DAOBeam]).timeFac]';
    timeFac(1) = [];
    timeFac(end) = [];
    %timeFac(timeFac == 0) = [];
    
    timeFacMatrix = sparse(i,j,timeFac,(apertureInfo.totalNumOfShapes-1),apertureInfo.totalNumOfShapes);
    timeBNOptAngles = timeFacMatrix*timeDAOBorderAngles;
    
    % values of average leaf speeds of optimized gantry angles
    c_lfspd = reshape([abs(diff(reshape(leftLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2)) ...
        abs(diff(reshape(rightLeafPos,apertureInfo.beam(1).numOfActiveLeafPairs,apertureInfo.totalNumOfShapes),1,2))]./ ...
        repmat(timeBNOptAngles',apertureInfo.beam(1).numOfActiveLeafPairs,2),2*apertureInfo.beam(1).numOfActiveLeafPairs*numel(timeBNOptAngles),1);
end

% values of doserate (MU/sec) in an arc surrounding the optimized angles
weights = apertureInfoVec(1:(apertureInfo.totalNumOfShapes))./apertureInfo.jacobiScale;
c_dosrt = apertureInfo.weightToMU.*weights./timeDoseBorderAngles;

% concatenate %Maybe do in a way we can call the superclass above?
c = [c_dos_dao; c_lfspd; c_dosrt];


