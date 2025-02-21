function stf = matRad_generateStf(ct,cst,pln,visMode)
% matRad steering information generation
%
% call
%   stf = matRad_generateStf(ct,cst,pln)
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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

generator = matRad_StfGeneratorBase.getGeneratorFromPln(pln);

if nargin == 4
    matRad_cfg.dispDeprecationWarning('The fourth ''visMode'' argument for matRad_generateStf is deprecated and will be removed in a future release. Please use pln.propStf.visMode as a replacement (or the corresponding property in the stf generators))');
    generator.visMode = visMode;
end

%call the calcDose funktion
stf = generator.generate(ct,cst);

end
