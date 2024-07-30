function matRad_visPhotonFieldShapes(pln)
% matRad function to visualize imported field shapes from a dicom RTPLAN
% data. Note that this only works with pln structures that were generated
% with matRad's dicom import tool and feature the appropriate field shape
% information
% 
% call
%   matRad_visPhotonFieldShapes(pln)
%
% input
%   pln: matRad plan struct
%
% output
%   -
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(pln, 'propStf') || ~isfield(pln.propStf, 'collimation')
    error(['This function only works with field shapes stemming from a DICOM ' ...
           'photon plan import. For matRad''s native aperture info struct you ' ...
           'may want to check out https://github.com/e0404/matRad/blob/master/matRad_visApertureInfo.m '...
           'as shown in https://github.com/e0404/matRad/blob/master/examples/matRad_example3_photonsDAO.m']);
end

coords = [-pln.propStf.collimation.fieldWidth/2 ...
         :pln.propStf.collimation.convResolution ...
         :pln.propStf.collimation.fieldWidth/2-pln.propStf.collimation.convResolution];

fieldShapeFigureHandle = figure; 
fieldShapeAxisHandle = axes(fieldShapeFigureHandle);

for i = 1: pln.propStf.collimation.numOfFields

    cla(fieldShapeAxisHandle)
    
    imagesc(fieldShapeAxisHandle, coords, coords, pln.propStf.collimation.Fields(i).Shape)
    
    xlabel(fieldShapeAxisHandle, '[mm]')
    ylabel(fieldShapeAxisHandle, '[mm]')

    title(fieldShapeAxisHandle, ['field shape i: ' num2str(i) ' @ gantry = ' num2str(pln.propStf.collimation.Fields(i).GantryAngle) '° and couch = ' num2str(pln.propStf.collimation.Fields(i).CouchAngle) '°'])
    axis(fieldShapeAxisHandle, 'equal', 'tight')
    
    drawnow
    pause(.1)
    
end