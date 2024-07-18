classdef matRad_None < matRad_BiologicalModel

    properties (Constant)
        model = 'none';
    end

    methods
        function this = matRad_None()
            this@matRad_BiologicalModel();
            this.availableRadiationModalities = {'photons', 'protons', 'carbon', 'helium'};
        end
    end
end