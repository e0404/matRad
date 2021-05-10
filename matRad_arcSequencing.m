function beam = matRad_arcSequencing(beam,stf,pln,weightToMU)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The sequencing algorithm generates an a priori unkown number of aperture.
% We only want to keep a certain number of them (numToKeep).  These will be
% the ones with the highest intensity-area product.
%
%
% call
%   beam =
%   matRad_arcSequencing(beam)
%
% input
%   beam:               beam struct with shapes and weights only at the
%                       initGantyAngles
%
%
% output
%   beam:               beam struct with shapes and weights distributed to
%                       the correct optGantryAngles
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


fileName = pln.propOpt.VMAToptions.machineConstraintFile;
try
    load(fileName,'machine');
catch
    error(['Could not find the following machine file: ' fileName ]);
end

numOfBeams = pln.propStf.numOfBeams;

leafDir = 1;

for i = 1:numOfBeams
    
    if stf(i).propVMAT.FMOBeam
        
        %Spread apertures to each child angle
        %according to the trajectory (mean leaf position). Assume that
        %shapes are already order in increased (left to right) position
        leafDir = -1*leafDir;
        
        childrenIndex = stf(i).propVMAT.beamChildrenIndex;
        if leafDir == -1
            % reverse order of shapes
            childrenIndex = flipud(childrenIndex);
        end
        
        count = 1;
        numOfShapes = beam(i).numOfShapes;
        
        for shape = 1:numOfShapes
            childIndex = childrenIndex(count);
            beam(childIndex).leafDir = leafDir;
            
            if childIndex == i
                % do not overwrite information, since we will need it for
                % the remaining beams (DAO, not init)
                beam(childIndex).tempNumOfShapes = 1;
                beam(childIndex).tempShapes = beam(i).shapes(:,:,shape);
                beam(childIndex).tempShapesWeight = beam(i).shapesWeight(shape);
                beam(childIndex).fluence = beam(childIndex).tempShapes;
                beam(childIndex).sum = beam(childIndex).tempShapesWeight*beam(childIndex).tempShapes;
            else
                % don't worry about overwriting
                beam(childIndex).numOfShapes = 1;
                beam(childIndex).shapes = beam(i).shapes(:,:,shape);
                beam(childIndex).shapesWeight = beam(i).shapesWeight(shape);
                beam(childIndex).fluence = beam(childIndex).shapes;
                beam(childIndex).sum = beam(childIndex).shapesWeight*beam(childIndex).shapes;
            end
            
            count = count+1;
        end
    else
        % if beam isn't an FMO beam, then there is no info in the beam
        % struct
        continue
    end
end

for i = 1:numOfBeams
    % now go through and calculate gantry rotation speed, MU rate, etc.
    
    if stf(i).propVMAT.FMOBeam
        beam(i).numOfShapes = beam(i).tempNumOfShapes;
        beam(i).shapes = beam(i).tempShapes;
        beam(i).shapesWeight = beam(i).tempShapesWeight;
        
        beam(i).tempNumOfShapes = [];
        beam(i).tempShapes = [];
        beam(i).tempShapesWeight = [];
        
        for j = 1:stf(i).propVMAT.numOfBeamSubChildren
            %Prevents matRad_sequencing2ApertureInfo from attempting to
            %convert shape to aperturevec for subchildren
            beam(stf(i).propVMAT.beamSubChildrenIndex(j)).numOfShapes = 0;
        end
    end
    
    if stf(i).propVMAT.DAOBeam
        beam(i).gantryRot = machine.constraints.gantryRotationSpeed(2); %gantry rotation rate until next opt angle
        beam(i).MURate = weightToMU.*beam(i).shapesWeight.*beam(i).gantryRot./stf(i).propVMAT.DAOAngleBordersDiff; %dose rate until next opt angle
        %Rescale weight to represent only this control point; weight will be shared
        %with the interpolared control points in matRad_daoVec2ApertureInfo
        beam(i).shapesWeight = beam(i).shapesWeight.*stf(i).propVMAT.timeFacCurr;
    end
end

beam = rmfield(beam,{'tempShapes','tempShapesWeight','tempNumOfShapes'});


