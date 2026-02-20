classdef matRad_EmptyBiologicalModel < matRad_BiologicalModel
%  matRad_None
%  Class to implement the None model
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Constant)
        model = 'none';
        possibleRadiationModes =  {'photons', 'protons', 'carbon', 'helium', 'brachy','VHEE'};
        requiredQuantities = {};
        defaultReportQuantity = 'physicalDose';
    end

    methods
        function this = matRad_EmptyBiologicalModel()
            this@matRad_BiologicalModel();
        end
    end
end