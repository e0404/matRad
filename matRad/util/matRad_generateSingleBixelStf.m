function stf = matRad_generateSingleBixelStf(ct,cst,pln)
% 
% call
%   stf = matRad_generateSingleBixelStf(ct,cst,pln,visMode)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   pln:        matRad plan meta information struct
%   visMode:    toggle on/off different visualizations by setting this value to 1,2,3 (optional)
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

matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispDeprecationWarning('The function %s is not intendet to be used, check out the Single Bixel Stf Generators in the steering folder',fileparts(mfilename));

if any(strcmp(pln.radiationMode,{'protons','helium','carbon','oxygen'}))
    stfGen = matRad_ParticleStfGeneratorSingleBeamlet(pln);
elseif strcmp(pln.radiationMode,'photons')
    stfGen = matRad_PhotonStfGeneratorSingleBeamlet(pln);
else
    matRad_cfg.dispError('Unsupported Radiation Mode!');
end

stf = stfGen.generate(ct,cst);

end    



