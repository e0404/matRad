function stf = matRad_generateStf(ct,cst,pln)
% matRad steering information generation
%
% call
%   stf = matRad_generateStf(ct,cst,pln,visMode)
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
    brachyStfGen = matRad_BrachyStfGenerator(pln);
    stf = brachyStfGen.generate(ct, cst);
elseif strcmp(pln.radiationMode, 'photons')
    photonStfGen = matRad_PhotonStfGenerator(pln);
    stf = photonStfGen.generate(ct, cst);
elseif any(strcmp(pln.radiationMode, {'protons', 'carbon', 'helium'}))
    ionStfGen = matRad_IonStfGenerator(pln);
    stf = ionStfGen.generate(ct, cst);
end

end