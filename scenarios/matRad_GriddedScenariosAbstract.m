classdef (Abstract) matRad_GriddedScenariosAbstract < matRad_ScenarioModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (AbortSet = true)
        includeNominalScenario = true;        
        combinations = 'none'; %Can be 'none', 'shift', 'all' to ontrol creation of worst case combinations 
        combineRange = true; %Wether to treat absolute & relative range as one shift or as separate scenarios
    end
    
    %Each subclass needs to define how many gridpoints it uses and if this
    %can be set or not
    properties (Abstract)
        numOfSetupGridPoints;
        numOfRangeGridPoints;
    end

    properties (Constant)
        validCombinationTypes = {'all','none','shift'};
    end

    methods
        function this = matRad_GriddedScenariosAbstract(ct)      
                         
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            %this.updateScenarios();
        end
        
        function set.combineRange(this,combineRange_)
            valid = isscalar(combineRange_) && (isnumeric(combineRange_) || islogical(combineRange_));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for combineRange! Needs to be a boolean / logical value!');
            end
            this.combineRange = combineRange_;
            this.updateScenarios();
        end

        %% set methods
        function set.includeNominalScenario(this,includeNomScen)
            valid = isscalar(includeNomScen) && (isnumeric(includeNomScen) || islogical(includeNomScen));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for includeNominalScenario! Needs to be a boolean / logical value!');
            end
            this.includeNominalScenario = includeNomScen;
            this.updateScenarios();
        end

        function set.combinations(this,combinations_)
            valid = any(strcmp(combinations_,this.validCombinationTypes));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for combinations! Needs to be one of the strings %s!',strjoin(this.validCombinationTypes,' / '));
            end
            this.combinations = combinations_;
            this.updateScenarios();
        end

        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();

            %Get the maximum, i.e., worst case shifts
            wcSetupShifts = this.wcSigma * this.shiftSD;
            
            %% Create gridded setup shifts
            %Create grid vectors for setup shifts
            setupShiftGrid = zeros(this.numOfSetupGridPoints,numel(wcSetupShifts));            
            if mod(this.numOfSetupGridPoints,2) == 0 && this.includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Setup Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end

            for i = 1:numel(wcSetupShifts)
                setupShiftGrid(:,i) = linspace(-wcSetupShifts(i),wcSetupShifts(i),this.numOfSetupGridPoints);
                if this.includeNominalScenario 
                      
                    [~,ix] = min(abs(setupShiftGrid(:,i)));
                    setupShiftGrid(ix,i) = 0;    
                end  
            end
            
            %Now create vector of all shifts for different combinatorial
            %settings
            switch this.combinations
                case 'none'
                    %Create independent shifts
                    griddedSetupShifts = [];
                    for i=1:size(setupShiftGrid,2)
                        tmpGrid = zeros(size(setupShiftGrid,1),3);
                        tmpGrid(:,i) = setupShiftGrid(:,i);
                        griddedSetupShifts = [griddedSetupShifts; tmpGrid];
                    end                    
                case {'shift','all'}
                    [X,Y,Z] = meshgrid(setupShiftGrid(:,1),setupShiftGrid(:,2),setupShiftGrid(:,3));
                    griddedSetupShifts = [X(:), Y(:), Z(:)];    
                otherwise
                    matRad_cfg.dispError('Invalid value for combinations! This sanity check should never be reached!');
            end

            griddedSetupShifts = matRad_ImportanceScenarios.uniqueStableRowsCompat(griddedSetupShifts);
            shiftNomScenIx = find(all(griddedSetupShifts == zeros(1,3),2));            
            
            if ~isempty(shiftNomScenIx) || this.includeNominalScenario
                if ~isempty(shiftNomScenIx)
                    griddedSetupShifts(shiftNomScenIx,:) = [];
                end
                griddedSetupShifts = [0 0 0; griddedSetupShifts];
            end
                        
            this.totNumShiftScen = size(griddedSetupShifts,1);
                                
            %% Create gridded range shifts
            %Obtain worst case range shifts
            wcRangeShifts = this.wcSigma * [this.rangeAbsSD this.rangeRelSD./100];        
            
            rangeShiftGrid = zeros(this.numOfRangeGridPoints,numel(wcRangeShifts));            
            if mod(this.numOfRangeGridPoints,2) == 0 && this.includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Range Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end

            for i = 1:numel(wcRangeShifts)
                rangeShiftGrid(:,i) = linspace(-wcRangeShifts(i),wcRangeShifts(i),this.numOfRangeGridPoints);
                if this.includeNominalScenario 
                      
                    [~,ix] = min(abs(rangeShiftGrid(:,i)));
                    rangeShiftGrid(ix,i) = 0;    
                end  
            end

            if this.combineRange
                griddedRangeShifts = rangeShiftGrid;
            else                
                [rngAbs,rngRel] = meshgrid(rangeShiftGrid(:,1),rangeShiftGrid(:,2));
                griddedRangeShifts = [rngAbs(:),rngRel(:)];
            end

            %Remove duplicate scenarios and update number of shifts
            griddedRangeShifts = this.uniqueStableRowsCompat(griddedRangeShifts); 

            rangeNomScenIx = find(all(griddedRangeShifts == zeros(1,2),2));            
            
            if ~isempty(rangeNomScenIx) || this.includeNominalScenario
                if ~isempty(rangeNomScenIx)
                    griddedRangeShifts(rangeNomScenIx,:) = [];
                end
                griddedRangeShifts = [0 0; griddedRangeShifts];
            end

            this.totNumRangeScen = size(griddedRangeShifts,1);
                       
            %Aggregate scenarios
            switch this.combinations
                case {'none','shift'}
                    scenarios = zeros(this.totNumShiftScen + this.totNumRangeScen,5);
                    scenarios(1:this.totNumShiftScen,1:3) = griddedSetupShifts;
                    scenarios(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,4:5) = griddedRangeShifts;

                    %create the linear mask of scenarios
                    linearMaskTmp = ones(size(scenarios,1),3);
                    linearMaskTmp(1:this.totNumShiftScen,2) = (1:this.totNumShiftScen)';
                    linearMaskTmp(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,3) = (1:this.totNumRangeScen)';

                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMaskTmp = linearMaskTmp(ia,:);

                case 'all'
                    %Prepare scenario matrix by replicating shifts
                    %with the number of range scenarios
                    scenarios = repmat(griddedSetupShifts,this.totNumRangeScen,1);
                    scenarios = [scenarios zeros(size(scenarios,1),2)];
                    
                    %create the linear mask of scenarios
                    linearMaskTmp = ones(size(scenarios,1),3);
                    for r = 1:this.totNumRangeScen                    
                        offset = (r-1)*this.totNumShiftScen;
                        ixR = (offset + 1):(offset + this.totNumShiftScen);
                        scenarios(ixR,4:5) = repmat(griddedRangeShifts(r,:),this.totNumShiftScen,1);

                        %Set linear mask
                        linearMaskTmp(ixR,2) = (1:this.totNumShiftScen)';
                        linearMaskTmp(ixR,3) = r;
                    end


                    %create the linear mask of scenarios
                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMaskTmp = linearMaskTmp(ia,:);
                otherwise
                    matRad_cfg.dispError('Invalid value for combinations! This sanity check should never be reached!');
            end

            %Handle 4D phases
            phases = repmat(1:this.numOfCtScen,size(scenarios,1),1);
            phases = phases(:);
            scenarios = horzcat(phases, repmat(scenarios,[this.numOfCtScen 1]));
            linearMaskTmp = repmat(linearMaskTmp,this.numOfCtScen,1);
            linearMaskTmp(:,1) = phases;

            %Finalize meta information
            this.totNumScen = size(scenarios,1);

            this.relRangeShift = scenarios(:,6);
            this.absRangeShift = scenarios(:,5);
            this.isoShift = scenarios(:,2:4);
            
            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

            this.scenMask = false(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);
            
            this.scenForProb = scenarios;
            this.linearMask = linearMaskTmp;
            
            maskIx = sub2ind(size(this.scenMask),linearMaskTmp(:,1),linearMaskTmp(:,2),linearMaskTmp(:,3));
            this.scenMask(maskIx) = true;

            %Get Scenario probability
            %First, we use the Gaussian Uncertainty model for range and
            %setup
            Sigma = diag([this.shiftSD,this.rangeAbsSD,this.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            
            tmpScenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios(:,2:end)/cs).^2, 2)) / prod(diag(cs));
            
            %Now we combine with the 4D ct phase probabilities (multiply)
            tmpPhaseProb = arrayfun(@(phase) this.phaseProbability(phase),phases);
            
            %Finalize probabilities
            this.scenProb = tmpPhaseProb .* tmpScenProb;
            this.scenWeight = this.scenProb./sum(this.scenProb);  

            %TODO: Discard scenarios with probability 0?
        end
    end

    methods (Static)
        function [uniqueStableRows,ia] = uniqueStableRowsCompat(values)
            %This is a compatability wrapper to call unique without sorting            
            
            matRad_cfg = MatRad_Config.instance();
            
            if matRad_cfg.isOctave
                %https://stackoverflow.com/questions/37828894/
                [~,ia,~] = unique(values,'rows','first');
                ia = sort(ia);
                uniqueStableRows = values(ia,:);
            else
                [uniqueStableRows,ia] = unique(values,'rows','stable');
            end
        end
    end
end