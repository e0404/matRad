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
        defaultRBEprotons = 1.1;
        defaultRBEphotons = 1;
        possibleRadiationModes = {'photons','protons','helium','carbon','brachy'};
        requiredQuantities = {};
    end

    properties
        RBE;
    end

    methods
        function this = matRad_ConstantRBE()
            this = this@matRad_BiologicalModel();

        end
  

        % function assignBioModelPropertiesFromEngine(this,engine)
        % 
        %     % This function mirrors the user defined property
        %     % pln.propDoseCalc.bioProperties.RBE to the model property RBE.
        %     % If not defined, just uses default values
        % 
        %     matRad_cfg = MatRad_Config.instance();
        % 
        %     if isprop(engine, 'bioProperties') && isfield(engine.bioProperties, 'RBE')
        %         this.RBE = engine.bioProperties.RBE;
        %     else
        %         switch radiationMode
        % 
        %             case 'photons'
        %                 this.RBE = this.defaultRBEphotons;
        %             case 'protons'
        %                 this.RBE = this.defaultRBEprotons;
        %         end
        % 
        %         matRad_cfg.dispWarning('No RBE value specified, using default value of %f', this.RBE);
        %     end
        % ends

        function calcAvailable = checkBioCalcConsistency(this, machine)

            matRad_cfg = MatRad_Config.instance();

            calcAvailable = checkBioCalcConsistency@matRad_BiologicalModel(this, machine);

            if isempty(this.RBE)
                matRad_cfg.dispWarning('No specifc constant RBE value provided, using default!');
                
                switch machine.meta.radiationMode
                    case 'photons'
                        this.RBE = this.defaultRBEphotons;

                    case 'protons'
                        this.RBE = this.defaultRBEprotons;
                end
            end
        end

    end


end