function sliceIndices = matRad_findRtssContourSlicesInCt(contourZ, ct)
% matRad function to find CT slices compatible with an RTSTRUCT contour plane
%
% call:
%   sliceIndices = matRad_findRtssContourSlicesInCt(contourZ,ct)
%
% input:
%   contourZ:       z-position of one RTSTRUCT contour plane in mm
%   ct:             matRad ct struct with z-axis and DICOM slice thickness
%
% output:
%   sliceIndices:   row vector with CT slice indices intersected by the
%                   physical slice slab of the contour plane
%
% References
%   -
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

matRad_cfg = MatRad_Config.instance();

if ~isnumeric(contourZ) || ~isscalar(contourZ) || ~isfinite(contourZ)
    matRad_cfg.dispError('contourZ must be a finite numeric scalar.');
end

sliceAxis = matRad_getCtSliceAxis(ct);
sliceThickness = matRad_getCtSliceThickness(ct);

halfSliceThickness = sliceThickness / 2;
tol = max(1e-6, 1e-6 * sliceThickness);

sliceIndices = find(contourZ + halfSliceThickness + tol > sliceAxis & ...
                    contourZ - halfSliceThickness - tol <= sliceAxis);
sliceIndices = sliceIndices(:)';

end

function sliceAxis = matRad_getCtSliceAxis(ct)

matRad_cfg = MatRad_Config.instance();

if isstruct(ct) && isfield(ct, 'z') && isnumeric(ct.z) && isvector(ct.z) && ...
        ~isempty(ct.z)
    sliceAxis = ct.z(:)';
elseif isstruct(ct) && isfield(ct, 'dicomInfo') && isstruct(ct.dicomInfo) && ...
        isfield(ct.dicomInfo, 'SlicePositions') && ...
        isnumeric(ct.dicomInfo.SlicePositions) && isvector(ct.dicomInfo.SlicePositions) && ...
        ~isempty(ct.dicomInfo.SlicePositions)
    sliceAxis = ct.dicomInfo.SlicePositions(:)';
else
    matRad_cfg.dispError('ct must contain a numeric z-axis or dicomInfo.SlicePositions.');
end

if any(~isfinite(sliceAxis))
    matRad_cfg.dispError('CT slice positions must be finite.');
end
end

function sliceThickness = matRad_getCtSliceThickness(ct)

matRad_cfg = MatRad_Config.instance();

if ~isstruct(ct) || ~isfield(ct, 'dicomInfo') || ~isstruct(ct.dicomInfo) || ...
        ~isfield(ct.dicomInfo, 'SliceThickness')
    matRad_cfg.dispError('ct.dicomInfo.SliceThickness is required.');
end

sliceThickness = ct.dicomInfo.SliceThickness;
if ~isnumeric(sliceThickness) || ~isscalar(sliceThickness) || ...
        ~isfinite(sliceThickness) || sliceThickness <= 0
    matRad_cfg.dispError('ct.dicomInfo.SliceThickness must be a positive numeric scalar.');
end
end
