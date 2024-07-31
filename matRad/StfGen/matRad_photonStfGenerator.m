classdef matRad_photonStfGenerator < matRad_externalStfGenerator

    properties (Constant)
        name = 'photonStfGen';
        shortName = 'photonStfGen';
    end 

    properties (Access = protected)
        
    end
    
    methods 
        function this = matRad_photonStfGenerator(pln)
            if nargin < 1
                pln = [];
            end
            this@matRad_externalStfGenerator(pln);
            this.radiationMode = 'photons';

            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));


            
            if ~isfield(pln, 'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end
         end
    end

    methods 

        function stf = generate(this, ct, cst)
            stf = generate@matRad_externalStfGenerator(this, ct, cst);
            this.initializePatientGeometry(ct, cst);
            stf = this.generateSourceGeometry(ct, cst);
        end
    end

    methods (Access = protected)
        function initializePatientGeometry(this, ct, cst)
            % Initialize the patient geometry
            initializePatientGeometry@matRad_externalStfGenerator(this);

        end
        
        function pbMargin = getPbMargin(this)
            pbMargin = getPbMargin@matRad_externalStfGenerator(this);
            pbMargin = pln.propStf.bixelWidth;
        end
        
        function stf = generateSourceGeometry(this, ct, cst)
            stf = generateSourceGeometry@matRad_externalStfGenerator(this);

            % book keeping for photons
            stf(i).ray(j).energy = machine.data.energy;
        end

        function initializeEnergy(this,ct,cst)
            initializeEnergy@matRad_externalStfGenerator(this, ct, cst);

            % book keeping for photons
            stf(i).ray(j).energy = machine.data.energy;
        end

    end
end
