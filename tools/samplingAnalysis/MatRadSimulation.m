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
        doseCubeDimensions
        nominalScenario
        ct
        cst
        pln
        
        % doseStat
        meanCube
        stdCube
        gammaCube
        gammaDoseAgreement
        gammaDistAgreement
        percentiles
        
        statisticsComputed = false
        
        % computed properties which need to be recalculated when dose
        % changes. is not in Dependent properties because calculation is
        % expensive
        dvhContainer
        qiContainer
    end
    
    properties (SetObservable, AbortSet, Access = private)
        % container of scenarios
        scenContainer
    end
    
    properties (Dependent = true)
        numOfScen
        weights
    end
    
    properties (Dependent = true, Access = private)
        nextScenIdentifier
        doseCubeDim
        weights_unnormalised
        subIx
    end
        
    methods
        function obj = MatRadSimulation(radiationQuantity, nominalScenario, ct, cst, pln, expectedNumOfScen)
            % optionals
            if ~exist('expectedNumOfScen', 'var') || isempty(expectedNumOfScen)
                expectedNumOfScen = 0;
            end
         
            obj.radiationQuantity = radiationQuantity;
            obj.nominalScenario = nominalScenario;
            obj.ct = ct;
            obj.cst = cst;
            obj.pln = pln;
            
            obj.scenContainer = NaN * ones(expectedNumOfScen, 1);
            obj.scenContainer = cell(1, expectedNumOfScen);
            addlistener(obj,'scenContainer','PostSet',@obj.handleChangeOfScen);
        end % eof constructor
        
        function obj = initNewScen(obj, scenario, w)
            if obj.numOfScen == 0
                obj.deriveConstantsFromFirstScenario(scenario);
            end
            if ~obj.isValidScen(scenario)
                error('Scenario not valid.');
            end
            
            % obj.weights_unnormalised(obj.nextScenIdentifier,1) = w;
            obj.scenContainer{obj.nextScenIdentifier} = scenario;
        end % eof initNewScen
        
        function obj = runAnalysis(obj, gammaCriteria, percentiles)
            obj.gammaDoseAgreement = gammaCriteria(1);
            obj.gammaDistAgreement = gammaCriteria(2);
            obj.percentiles = percentiles;
            
            obj.fillDoseContainer();
            obj.computeMeanStdCube();
            obj.computeGammaCube();   
            obj.computeAllDvhQi();
            obj.computeStructureWiseDvhQi();
            
        end
        
        function weights_unnormalised = get.weights_unnormalised(obj)
          weights_unnormalised = NaN * ones(obj.numOfScen, 1);
            for i = 1:obj.numOfScen
                weights_unnormalised(i) = obj.scenContainer{i}.weight;
            end
        end % eof get.NumOfScen
        
        function weights = get.weights(obj)
            weights = obj.weights_unnormalised ./ sum(obj.weights_unnormalised);
        end % eof get.NumOfScen
        
        function numOfScen = get.numOfScen(obj)
            numOfScen = sum(~cellfun(@isempty,obj.scenContainer));
        end % eof get.NumOfScen
        
        function doseCubeDim = get.doseCubeDim(obj)
            doseCubeDim = obj.nominalScenario.doseCubeDim;
        end % eof get.doseCubeDim
    
        function subIx = get.subIx(obj)
            subIx = obj.nominalScenario.subIx;
        end % eof get.subIx
        
        function numOfScen = get.nextScenIdentifier(obj)
            numOfScen = obj.numOfScen + 1;
        end % eof get.NumOfScen
    end
    
    methods (Static)
        
    end

    methods (Access = private)
         function computeAllDvhQi(obj, refVol, refGy)
            if ~exist('refVol', 'var') || isempty(refVol) || ~exist('refGy', 'var') || isempty(refVol)
              refVol = [2 5 50 95 98];
              refGy = linspace(0,max(obj.nominalScenario.dose(:)),6);
            end
            % compute nomQIDVH
            obj.nominalScenario.calcQiDVH(obj.cst, obj.pln, 'cum', [], refGy, refVol);
            doseGrid = obj.nominalScenario.dvh(1).doseGrid;
            for i = 1:obj.numOfScen
                obj.scenContainer{i}.calcQiDVH(obj.cst, obj.pln, 'cum', doseGrid, refGy, refVol);
            end
        end
        
        function computeStructureWiseDvhQi(obj)
          dvhContainer = struct();
          qiContainer = struct();
          % reassing dvh to stats structure
          for i = 1:size(obj.cst,1)
              dvhContainer(i).name = obj.cst{i,2};
              qiContainer(i).name = obj.cst{i,2};
              dvhContainer(i).weights         = obj.weights;
              qiContainer(i).weights          = obj.weights;
              dvh = obj.scenContainer{1}.getSingleStructDVH(dvhContainer(i).name);
              dvhContainer(i).doseGrid        = NaN * ones(obj.numOfScen, numel(dvh.doseGrid));
              dvhContainer(i).volumePoints    = NaN * ones(obj.numOfScen, numel(dvh.volumePoints));
              for j = 1:obj.numOfScen
                  dvh = obj.scenContainer{j}.getSingleStructDVH(dvhContainer(i).name);
                  dvhContainer(i).doseGrid(j,:)     = dvh.doseGrid;
                  dvhContainer(i).volumePoints(j,:) = dvh.volumePoints;
                  qiContainer(i).qi(j)                = obj.scenContainer{j}.getSingleStructQi(qiContainer(i).name);
              end  
          end
          obj.dvhContainer = dvhContainer;
          obj.qiContainer = qiContainer;
        end
        
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
            obj.dvhContainer = struct();
            obj.qiContainer = struct();
            obj.doseContainer = [];
            
            fprintf('Changed.');
        end
        
        function obj = fillDoseContainer(obj)
            for i=1:obj.numOfScen
                obj.doseContainer(:,i) = obj.scenContainer{i}.dose(:);
            end
        end
        
        function computeMeanStdCube(obj)
          obj.meanCube = zeros(obj.doseCubeDim);
          obj.stdCube  = zeros(obj.doseCubeDim);
          
          ix = obj.subIx;
          
          obj.meanCube(ix) = (sum(obj.doseContainer * diag(obj.weights),2));
          obj.stdCube(ix)  = std(obj.doseContainer, obj.weights,2);
        end
        
        function computeGammaCube(obj)
          obj.gammaCube = zeros(obj.doseCubeDim);

          obj.gammaCube = matRad_gammaIndex(obj.nominalScenario.dose,obj.meanCube,[obj.ct.resolution.x obj.ct.resolution.y obj.ct.resolution.z],[obj.gammaDoseAgreement obj.gammaDistAgreement]);
        end
              
    end
    
end
