classdef matRad_ConstantRBE < matRad_BiologicalModel
%  matRad_ConstantRBE
%  Class to implement the constantRBE model
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
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        model = 'constRBE'
        possibleRadiationModes = {'photons','protons','helium','carbon','brachy'};
        requiredQuantities = {'physicalDose'};
        defaultReportQuantity = 'RBExDose';
    end

    properties
        RBE = 1.1;
    end

    methods
        function this = matRad_ConstantRBE()
            this = this@matRad_BiologicalModel();
        end
    end


end