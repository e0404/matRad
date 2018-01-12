classdef MatRadSimulation < handle
    
    properties (Access = private)
        % computed properties which need to be recalculated when dose
        % changes. is not in Dependent properties because calculation is
        % expensive
        doseContainer
    end

    properties (SetAccess = private)
        % should only be changed with constructor
        radiationQuantity
        subIx
        doseCubeDimensions
        nominalScenario
        cst
        
        statisticsComputed = false
    end
    
    properties (SetObservable, AbortSet, Access = private)
        % container of scenarios
        scenContainer
    end
    
    properties (Dependent = true)
        numOfScen
        dvhContainer
        qiContainer
    end
    
    properties (Dependent = true, Access = private)
        nextScenIdentifier
    end
        
    methods
        function obj = MatRadSimulation(radiationQuantity, nominalScenario, cst, subIx, expectedNumOfScen)
            % optionals
            if ~exist('expectedNumOfScen', 'var') || isempty(expectedNumOfScen)
                expectedNumOfScen = 0;
            end
         
            obj.radiationQuantity = radiationQuantity;
            obj.nominalScenario = nominalScenario;
            obj.nominalScenario = cst;
            obj.subIx = subIx;
            
            obj.scenContainer = cell(1, expectedNumOfScen);
            addlistener(obj,'scenContainer','PostSet',@obj.handleChangeOfScen);
        end % eof constructor
        
        function obj = initNewScen(obj, scenario)
            if obj.numOfScen == 0
                obj.deriveConstantsFromFirstScenario(scenario);
            end
            if ~obj.isValidScen(scenario)
                error('Scenario not valid.');
            end
            
            obj.scenContainer{obj.nextScenIdentifier} = scenario;
        end % eof initNewScen
        
        function obj = runAnalysis(obj)
            obj.fillDoseContainer();
        end
        
        function numOfScen = get.numOfScen(obj)
            numOfScen = sum(~cellfun(@isempty,obj.scenContainer));
        end % eof get.NumOfScen
    
        function numOfScen = get.nextScenIdentifier(obj)
            numOfScen = obj.numOfScen + 1;
        end % eof get.NumOfScen
    end
    
    methods (Static)
        
    end

    methods (Access = private)
        function validScen = isValidScen(obj, scenario)
            if scenario.radiationQuantity ~= obj.radiationQuantity
                error('Scenarios can only be added if they are of the same quantity.');
            else
                validScen = true;
            end
            
        end % eof isValidScen
        
        function obj = deriveConstantsFromFirstScenario(obj, scenario)
            obj.doseCubeDimensions = size(scenario.dose);            
        end % eof deriveConstantsFrom
        
        function obj = handleChangeOfScen(obj,~,~)
            obj.statisticsComputed = false;
            fprintf('Changed.');
        end
        
        function obj = fillDoseContainer(obj)
            for i=1:obj.numOfScen
                obj.doseContainer(:,i) = obj.scenContainer{i}.dose(:);
            end
        end
              
    end
    
end
