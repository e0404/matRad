classdef matRad_NominalScenario < matRad_ScenarioModel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        name = 'nomScen';
    end

    methods
        function this = matRad_NominalScenario(ct)
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end            
            this@matRad_ScenarioModel(superclassArgs{:});
        end
        
        function scenarios = updateScenarios(this)
            
            %Scenario Probability from pdf - here it is one since only one
            %scenario exist
            %TODO: In the context of an uncertainty model, we should
            %consider assigning probability according to the model, and
            %just leaving the weight 1
            this.scenForProb = [0 0 0 0 0];
            this.scenProb = 1;

            %Scenario weight 
            this.scenWeight = 1;

            %set variables
            this.totNumShiftScen = 1;
            this.totNumRangeScen = 1;
            this.totNumScen = this.numOfCtScen; 
            
            %Individual shifts
            this.relRangeShift = 0;
            this.absRangeShift = 0;
            this.isoShift = [0 0 0];

            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.absRangeShift);

            %Mask for scenario selection
            this.scenMask = true(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);
            
            %generic code
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            totNumScen    = sum(this.scenMask(:));

            if totNumScen ~= this.totNumScen
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
        
    end
end

