function multScen = matRad_multScen(ct,scenarioModel)
%  matRad_multScen
%  This function replaces the deprecated matRad_multScen class. It creates
%  the respective matRad_ScenarioModel instance for the specific scenario
%  model chosen with standard parameters. 
%
% call
%   pln.multScen = matRad_multScen(ct,scenarioModel);
%
%   e.g. pln.multScen = matRad_multScen(ct,'nomScen');
%
% input
%   ct:                 ct cube
%   scenarioModel:      string to denote a scenario creation method
%                       'nomScen'   create only the nominal scenario
%                       'wcScen'    create worst case scenarios
%                       'impScen'   create important/grid scenarios
%                       'rndScen'   create random scenarios
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matRad_cfg = MatRad_Config.instance();
%matRad_cfg.dispWarning('The matRad_multScen function will be deprecated soon!\nCheck out the new Scenario Models in the scenarios folder.');

multScen = matRad_ScenarioModel.create(scenarioModel,ct);

end









