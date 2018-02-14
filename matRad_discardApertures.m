function beam = matRad_discardApertures(beam,numToKeep)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The sequencing algorithm generates an a priori unkown number of aperture.
% We only want to keep a certain number of them (numToKeep).  These will be
% the ones with the highest intensity-area product.
%
%
% call
%   beam =
%   matRad_discardApertures(beam,numToKeep)
%
% input
%   beam:               beam struct containing original shapes and
%                       intensities
%
%   numToKeep:          number of apertures to keep
%
% output
%   beam:               beam struct containing shapes and re-scaled
%                       intensities of the apertures we are keeping
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


%Find the numToKeep apertures having the highest dose-area product
numToKeep = min(numToKeep,k);
for shape = 1:k
    beam.DAP(shape) = (stf.bixelWidth)^2.*nnz(beam.shapes(:,:,shape)).*beam.shapesWeight(shape);
    x = repmat(1:size(beam.shapes(:,:,shape),2),size(beam.shapes(:,:,shape),1),1);
    comPosRow = sum(beam.shapes(:,:,shape).*x,2)./sum(beam.shapes(:,:,shape),2);
    beam.comPos(shape) = mean(comPosRow(~isnan(comPosRow),1));
end

[~,beam.comPosToDAPSort] = sort(beam.DAP,'descend');

totDAP_all = sum(beam.DAP(:));
totDAP_keep = sum(beam.DAP(beam.comPosToDAPSort(1:numToKeep)));


tempShapes = zeros(dimOfFluenceMxZ,dimOfFluenceMxX,numToKeep);
tempShapesWeight = zeros(numToKeep,1);
tempNewDAP = tempShapesWeight;
tempComPos = tempShapesWeight;
segmentKeep = 1;

%Keep only those numToKeep apertures with the highest DAP
%Preserve the shapes of the apertures, but scale the weights so
%that the total DAP is kept
for shape = 1:k
    if beam.comPosToDAPSort(shape) <= numToKeep
        tempShapes(:,:,segmentKeep) = beam.shapes(:,:,shape);
        tempNewDAP(segmentKeep) = totDAP_all*beam.DAP(shape)/totDAP_keep;
        tempShapesWeight(segmentKeep) = tempNewDAP(segmentKeep)/((stf.bixelWidth)^2.*nnz(tempShapes(:,:,segmentKeep))); %sequencing.beam.shapesWeight(sequencing.beam.segmentSortedDAP(segment))
        tempComPos(segmentKeep) = beam.comPos(shape);
        
        segmentKeep = segmentKeep+1;
    else
        continue
    end
end
beam.numOfShapes  = numToKeep;
beam.shapes = tempShapes;
beam.shapesWeight = tempShapesWeight;
beam.DAP = tempNewDAP;
beam.comPos = tempComPos;


