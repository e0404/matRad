classdef matRad_ConstantRBE < matRad_BiologicalModel

    properties (Constant)

        model = 'constRBE'
        defaultRBEprotons = 1.1;
        defaultRBEphotons = 1;
        

    end

    properties (SetAccess = protected)
        RBE;
    end

    methods
        function this = matRad_ConstantRBE()
            this = this@matRad_BiologicalModel();
            this.availableRadiationModalities = {'protons', 'photons'};
        end
  

        function assignBioModelPropertiesFromPln(this,engine)
            
            matRad_cfg = MatRad_Config.instance();

            if isprop(engine, 'bioProperties') && isfield(engine.bioProperties, 'RBE')
                this.RBE = engine.bioProperties.RBE;
            else
                switch radiationMode
  
                    case 'photons'
                        this.RBE = this.defaultRBEphotons;
                    case 'protons'
                        this.RBE = this.defaultRBEprotons;
                end
                
                matRad_cfg.dispWarning('No RBE value specified, using default value of %f', this.RBE);
            end
        end

    end
end