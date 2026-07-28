function [xLeft, xRight, yRow] = matRad_apertureQuads(apertureInfo, i, leftLeafPos, rightLeafPos, leafCoord, bixelWidth)
% matRad_apertureQuads Vertices of the closed parts of the leaf pairs.
%
%   Returned as 4-by-nLeafPairs matrices - one column per quadrangle, i.e.
%   one patch/fill face per leaf. Taking the leaf positions as an argument
%   rather than reading them from a shape lets the aperture player draw
%   positions interpolated between two shapes.
%
% call:
%   [xLeft, xRight, yRow] = matRad_apertureQuads(apertureInfo, i, ...
%                                                leftLeafPos, rightLeafPos, ...
%                                                leafCoord, bixelWidth)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minX = apertureInfo.beam(i).MLCWindow(1);
maxX = apertureInfo.beam(i).MLCWindow(2);
nK = apertureInfo.beam(i).numOfActiveLeafPairs;

[rowLo, rowHi] = matRad_leafRowBounds(apertureInfo, i, leafCoord, bixelWidth);
yRow = [rowLo; rowLo; rowHi; rowHi];

leftLeafPos = reshape(leftLeafPos(1:nK), 1, []);
rightLeafPos = reshape(rightLeafPos(1:nK), 1, []);

xLeft = [repmat(minX, 1, nK); leftLeafPos; leftLeafPos; repmat(minX, 1, nK)];
xRight = [rightLeafPos; repmat(maxX, 1, nK); repmat(maxX, 1, nK); rightLeafPos];

end

function [rowLo, rowHi] = matRad_leafRowBounds(apertureInfo, i, leafCoord, bixelWidth)
% Lower/upper edge of every active leaf row of beam i, as row vectors.

nK = apertureInfo.beam(i).numOfActiveLeafPairs;

if strcmp(leafCoord, 'physical')
    rowCentre = apertureInfo.beam(i).leafPairPos(1:nK);
    halfRow = bixelWidth / 2;
else
    % the leaf indices have to be flipped in order to fit to the order of
    % the leaf positions (1st row of leafPos is lowest row in physical
    % coordinates)
    rowCentre = flipud(find(apertureInfo.beam(i).isActiveLeafPair));
    rowCentre = rowCentre(1:nK);
    halfRow = 1 / 2;
end

rowCentre = reshape(rowCentre, 1, []);
rowLo = rowCentre - halfRow;
rowHi = rowCentre + halfRow;

end
