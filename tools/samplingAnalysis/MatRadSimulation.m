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
        meanDoseScenario
        stdDoseScenario
        gammaDoseScenario
        
        gammaDoseAgreement
        gammaDistAgreement
        percentiles
        
        statisticsComputed = false
        
        % computed properties which need to be recalculated when dose
        % changes. is not in Dependent properties because calculation is
        % expensive
        dvhDoseGrid
        dvhContainer
        qiContainer
        
        dvhStatistics
        qiStatistics
    end
    
    properties (SetObservable, AbortSet, SetAccess = private)
        % container of scenarios
        scenContainer
    end
    
    properties (Dependent = true)
        numOfScen
        weights
        metric
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
            obj.scenContainer = cell(expectedNumOfScen, 1);
            addlistener(obj,'scenContainer','PostSet',@obj.handleChangeOfScen);
        end % eof constructor
        
        function obj = initNewScen(obj, scenario)
            if obj.numOfScen == 0
                obj.deriveConstantsFromFirstScenario(scenario);
            end
            if ~obj.isValidScen(scenario)
                error('Scenario not valid.');
            end
            
            % obj.weights_unnormalised(obj.nextScenIdentifier,1) = w;
            obj.scenContainer{obj.nextScenIdentifier} = scenario;
        end % eof initNewScen
        
        function obj = initFractionatedTreatments(obj, scenarios, treatmentNumber, samplingMethod)
            if ~strcmp(samplingMethod, 'random')
                error('Only random sampling possible at the moment.')
            end
            if obj.numOfScen == 0
                obj.deriveConstantsFromFirstScenario(scenarios{1});
            end
            if size(treatmentNumber,1) ~= 1 % must be row vector to be iterable
                treatmentNumber = treatmentNumber';
            end
            treatmentNumber_unique = unique(treatmentNumber);

            % go through complete treatment scenarios
            for treatment_currIx = treatmentNumber_unique
                % select all scenarios which belong to current treatment_ix
                scenForCurrentTreatment = find(treatmentNumber == treatment_currIx);
                cumDose = zeros('like',scenarios{1}.doseLin);
                ctIndex = NaN * ones(numel(scenForCurrentTreatment),1);
                shifts = NaN * ones(numel(scenForCurrentTreatment),3);
                relRangeShifts =  NaN * ones(numel(scenForCurrentTreatment),1);
                absRangeShifts = NaN * ones(numel(scenForCurrentTreatment),1);
                runIx = 1;
                for s = scenForCurrentTreatment
                    if ~obj.isValidScen(scenarios{s})
                        error('Scenario not valid.');
                    end
                    cumDose = cumDose + scenarios{s}.doseLin;
                    ctIndex(runIx) = scenarios{s}.ctShiftIdentifier;
                    shifts(runIx, :) = scenarios{s}.shift;
                    relRangeShifts(runIx) = scenarios{s}.relRangeShift;
                    absRangeShifts(runIx) = scenarios{s}.absRangeShift;
                    runIx = runIx + 1;
                end
                obj.scenContainer{obj.nextScenIdentifier} = MatRadScenario(cumDose, scenarios{1}.subIx, 'phsicalDose', ctIndex, shifts, ...
                            relRangeShifts, absRangeShifts, scenarios{1}.doseCubeDim, 1);
            end
        end
        
        function obj = runAnalysis(obj, gammaCriteria, percentiles)
            obj.gammaDoseAgreement = gammaCriteria(1);
            obj.gammaDistAgreement = gammaCriteria(2);
            obj.percentiles = percentiles;
            
            obj.fillDoseContainer();
            obj.computeMeanStdCube();
            obj.computeGammaCube();   
            obj.computeAllDvhQi();
            obj.computeStructureWiseDvhQi();
            obj.computeStructureWiseStatistics();
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
        
        function metric = get.metric(obj)
            percentileNames = cell(numel(obj.percentiles),1);
            % create fieldnames
            for i = 1:numel(obj.percentiles)
                percentileNames{i} = ['P',num2str(obj.percentiles(i)*100)];
            end
            % create table rownames
            metric = vertcat({'mean';'min';'max';'std'},percentileNames{:});
        end
    end
    
    methods (Static)
        function S = wMean(X,w)
          if exist('w','var') || ~isempty(w)
              if isvector(X) && isvector(w)
                  S = reshape(w,1,[]) * reshape(X,[],1) / sum(w);
              else
                  % row-wise
                  S = reshape(w,1,[]) * X ./ sum(w);        
              end

          else
              S = mean(X);
          end
        end % eof wMean
        

        

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
            obj.dvhDoseGrid = doseGrid;
        end
        
        function computeStructureWiseDvhQi(obj)
          obj.dvhContainer = struct();
          obj.qiContainer = struct();
          % reassing dvh to stats structure
          for i = 1:size(obj.cst,1)
              obj.dvhContainer(i).name = obj.cst{i,2};
              obj.qiContainer(i).name = obj.cst{i,2};
              %obj.dvhContainer(i).weights         = obj.weights;
              %obj.qiContainer(i).weights          = obj.weights;
              dvh = obj.scenContainer{1}.getSingleStructDVH(obj.dvhContainer(i).name);
              obj.dvhContainer(i).doseGrid        = NaN * ones(1, numel(dvh.doseGrid));
              obj.dvhContainer(i).volumePoints    = NaN * ones(obj.numOfScen, numel(dvh.volumePoints));
              for j = 1:obj.numOfScen
                  dvh = obj.scenContainer{j}.getSingleStructDVH(obj.dvhContainer(i).name);
                  if j == 1
                      obj.dvhContainer(i).doseGrid(1,:) = dvh.doseGrid;
                  else
                      if obj.dvhContainer(i).doseGrid(1,:) ~= dvh.doseGrid
                          error('Dose grids are not equal.')
                      end
                  end
                  obj.dvhContainer(i).volumePoints(j,:) = dvh.volumePoints;
                  obj.qiContainer(i).qi(j)                = obj.scenContainer{j}.getSingleStructQi(obj.qiContainer(i).name);
              end  
          end
        end
        
        function computeStructureWiseStatistics(obj)
            % create statstics where structure based results (QI and DVH) are available
            for i = 1:numel(obj.dvhContainer)
              obj.dvhStatistics(i).VOIname     = obj.dvhContainer(i).name;
                obj.dvhStatistics(i).stat     = obj.calcDVHStat(obj.dvhContainer(i).volumePoints, obj.percentiles, obj.weights);
              
              obj.qiStatistics(i).name      = obj.qiContainer(i).name;
              obj.qiStatistics(i).stat      = obj.calcQiStat(obj.qiContainer(i).qi, obj.percentiles, obj.weights);
            end
        end
        
        function validScen = isValidScen(obj, scenario)
            if ~strcmp(scenario.radiationQuantity,obj.radiationQuantity)
                error('Scenarios can only be added if they are of the same quantity.');
                validScen = false;
            else
                validScen = true;
            end
            
            if ~(scenario.subIx == obj.subIx)
                error('Scenarios can only be added if they feature the same sub indices');
                validScen = false;
            else
                validScen = true;
            end
            
            if ~(scenario.doseCubeDim == obj.doseCubeDimensions)
                error('Scenarios can only be added if their cube dimension agree.');
                validScen = false;
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
            obj.dvhStatistics = struct();
            obj.qiContainer = struct();
            obj.qiStatistics = struct();
            obj.doseContainer = [];
            
            fprintf('Changed.');
        end
        
        function obj = fillDoseContainer(obj)
            for i=1:obj.numOfScen
                obj.doseContainer(:,i) = obj.scenContainer{i}.doseLin;
            end
        end
        
        function computeMeanStdCube(obj)
          meanCube = zeros(obj.doseCubeDim);
          stdCube  = zeros(obj.doseCubeDim);
          
          ix = obj.subIx;
          
          meanCube(ix) = (sum(obj.doseContainer * diag(obj.weights),2));
          stdCube(ix)  = std(obj.doseContainer, obj.weights,2);
          
          obj.meanDoseScenario = MatRadScenario(meanCube, obj.subIx, obj.radiationQuantity, NaN, NaN, NaN, NaN, obj.doseCubeDim, NaN)
          obj.stdDoseScenario = MatRadScenario(stdCube, obj.subIx, obj.radiationQuantity, NaN, NaN, NaN, NaN, obj.doseCubeDim, NaN)
        end
        
        function computeGammaCube(obj)
          gammaCube = zeros(obj.doseCubeDim);
          gammaCube = matRad_gammaIndex(obj.nominalScenario.dose,obj.meanDoseScenario.dose,[obj.ct.resolution.x obj.ct.resolution.y obj.ct.resolution.z],[obj.gammaDoseAgreement obj.gammaDistAgreement]);
          
          obj.gammaDoseScenario = MatRadScenario(gammaCube, obj.subIx, obj.radiationQuantity, NaN, NaN, NaN, NaN, obj.doseCubeDim, NaN);
        end
        
        function dvhStat = calcDVHStat(obj, volumePoints, percentiles,w)
            dvhMat = volumePoints; % rows are different scenarios, columns positions
            % for statistical reasons, treat NaN as 0
            dvhMat(isnan(dvhMat)) = 0;
            
            dvhStat.mean.volumePoints = obj.wMean(dvhMat,w);
            dvhStat.min.volumePoints = min(dvhMat);
            dvhStat.max.volumePoints = max(dvhMat);
            dvhStat.std.volumePoints = std(dvhMat,w);
    
            dvhStat.percDVH = NaN * ones(numel(percentiles),size(volumePoints, 2));
            
            for j = 1:size(dvhMat,2)
                wQ =  matRad_weightedQuantile(dvhMat(:,j), percentiles, w', false, 'none');
                dvhStat.percDVH(:,j) = wQ;
            end
    
        end % eof calcDVHStat
        
        function qiStat = calcQiStat(obj, qi,percentiles,w)
            fields = fieldnames(qi);
            % remove name field
            if sum(strcmp('name', fields)) >= 1
                qi = rmfield(qi, 'name');
            end
            fields = fieldnames(qi);
            qiStruct = qi;
            
            % create helper matlab structure which will be converted to table
            qiStatH = struct();
            for j = 1:numel(fields)
                if numel([qiStruct(:).(fields{j})]) == numel(w)
                    qiStatH(1).(fields{j}) = obj.wMean([qiStruct(:).(fields{j})],w);
                    qiStatH(2).(fields{j}) = min([qiStruct(:).(fields{j})]);
                    qiStatH(3).(fields{j}) = max([qiStruct(:).(fields{j})]);
                    qiStatH(4).(fields{j}) = std([qiStruct(:).(fields{j})],w);
                    wQ = matRad_weightedQuantile([qiStruct(:).(fields{j})], percentiles, w', false, 'none');
                    for k = 1:numel(wQ)
                        sIx = k + 4;
                        qiStatH(sIx).(fields{j}) = wQ(k);
                    end
                else
                    for k = 1:(4 + numel(percentiles))
                        qiStatH(k).(fields{j}) = [];
                    end
                end
            end
            qiStat = struct2table(qiStatH);
            qiStat.Properties.RowNames = obj.metric;
        end % eof calcQiStat
    end
    
end
