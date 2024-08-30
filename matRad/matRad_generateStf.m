% matRad steering information generation
%
%
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting
%               this value to 1,2,3 (optional)
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(pln.radiationMode, 'brachy')
    brachyStfGen = matRad_brachyStfGenerator(pln);
    stf = brachyStfGen.generate(ct, cst, 1);
elseif strcmp(pln.radiationMode, 'photons')
    photonStfGen = matRad_photonStfGenerator(pln);
    stf = photonStfGen.generate(ct, cst, 1);
elseif any(strcmp(pln.radiationMode, {'protons', 'carbon', 'helium'}))
    ionStfGen = matRad_ionStfGenerator(pln);
    stf = ionStfGen.generate(ct, cst, 1);
end
