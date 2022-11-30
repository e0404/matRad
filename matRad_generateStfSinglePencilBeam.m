function stf = matRad_generateStfSinglePencilBeam(ct,cst,pln)
% matRad steering information generation
% 
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%
% output
%   stf:        matRad steering information struct
%
% References
%   -
%
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


matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('matRad: Generating a single pencil beam stf struct...\n');

if (numel(pln.propStf.gantryAngles) ~= numel(pln.propStf.couchAngles)) && (numel(pln.propStf.gantryAngles) > 1)
    matRad_cfg.dispError('Inconsistent number of gantry and couch angles.');
end

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep 'basedata' filesep fileName]);
   SAD = machine.meta.SAD;
catch
   matRad_cfg.dispError('Could not find the following machine file: %s',fileName); 
end
    
% Define steering file like struct. Prellocating for speed.
stf = struct;

% Save meta information for treatment plan
stf.gantryAngle   = pln.propStf.gantryAngles;
stf.couchAngle    = pln.propStf.couchAngles;
stf.bixelWidth    = pln.propStf.bixelWidth;
stf.radiationMode = pln.radiationMode;
stf.SAD           = SAD;
stf.isoCenter     = pln.propStf.isoCenter(1,:);
stf.numOfRays     = 1;

% Save ray and target position in beam eye's view (bev)
stf.ray.rayPos_bev = [0 0 0];
stf.ray.targetPoint_bev = [0 SAD 0];

% source position in bev
stf.sourcePoint_bev = [0 -SAD 0];

% get (active) rotation matrix 
% transpose matrix because we are working with row vectors
rotMat_vectors_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles,pln.propStf.couchAngles));


stf.sourcePoint = stf.sourcePoint_bev*rotMat_vectors_T;

% Save ray and target position in lps system.  
stf.ray.rayPos      = stf.ray.rayPos_bev*rotMat_vectors_T;
stf.ray.targetPoint = stf.ray.targetPoint_bev*rotMat_vectors_T;
   

% loop over all rays to determine meta information for each ray    
stf.numOfBixelsPerRay = 1;

stf.ray.energy = 100;

stf.ray.rangeShifter.ID = 0;
stf.ray.rangeShifter.eqThickness = 0;
stf.ray.rangeShifter.sourceRashiDistance = 0;

stf.ray.focusIx = 1;  
          
          
% save total number of bixels
stf.totalNumOfBixels = 1;

end
