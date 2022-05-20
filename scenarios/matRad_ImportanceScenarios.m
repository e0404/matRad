classdef matRad_ImportanceScenarios < matRad_ScenarioModel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        includeNominalScenario = true;        

        numOfSetupGridPoints = 9;
        numOfRangeGridPoints = 9;
    end

    properties (SetAccess=protected)
        name = 'impScen';
    end

    properties (Access = private)
        combinations = 'none'; %Can be 'none', 'shift', 'all' to ontrol creation of worst case combinations 
        combineRange = true; %Wether to treat absolute & relative range as one shift or as separate scenarios
    end
    
    methods
        function this = matRad_ImportanceScenarios(ct)           
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});
        end
        
        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();
              
            %Get the maximum, i.e., worst case shifts
            wcSetupShifts = this.wcSigma * this.shiftSD;
            
            %% Create gridded setup shifts
            %Create grid vectors for setup shifts
            setupShiftGrid = zeros(this.numOfSetupGridPoints,numel(wcSetupShifts));            
            if mod(this.numOfSetupGridPoints,2) == 0
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
                case 'shift'
                    error('Not implemented!');
                    %Remove duplicate scenarios and update number of shifts
                    %griddedSetupShifts = unique(griddedSetupShifts,'rows','stable');    
                case 'all'
                    error('Not implemented!');
            end

            griddedSetupShifts = unique(griddedSetupShifts,'rows','stable');
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
            if mod(this.numOfRangeGridPoints,2) == 0
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
                matRad_cfg.dispError('Range Shift Permutations are not yet implemented!');
                
                %Preliminary code
                for i=1:size(setupShiftGrid,2)
                    tmpGrid = zeros(size(setupShiftGrid,1),3);
                    tmpGrid(:,i) = setupShiftGrid(:,i);
                    griddedRangeShifts = [griddedRangeShifts; tmpGrid];
                end
            end

            %Remove duplicate scenarios and update number of shifts
            griddedRangeShifts = unique(griddedRangeShifts,'rows','stable'); 

            rangeNomScenIx = find(all(griddedRangeShifts == zeros(1,2),2));            
            
            if ~isempty(rangeNomScenIx) || this.includeNominalScenario
                if ~isempty(rangeNomScenIx)
                    griddedRangeShifts(rangeNomScenIx,:) = [];
                end
                griddedRangeShifts = [0 0; griddedRangeShifts];
            end

            this.totNumRangeScen = size(griddedRangeShifts,1);
            
            
            %Aggregate scenarios

            if strcmp(this.combinations,'none')
                scenarios = zeros(this.totNumShiftScen + this.totNumRangeScen,5);
                scenarios(1:this.totNumShiftScen,1:3) = griddedSetupShifts;
                scenarios(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,4:5) = griddedRangeShifts;
                

                
                %create the linear mask of scenarios
                linearMask = ones(size(scenarios,1),3);
                linearMask(1:this.totNumShiftScen,2) = (1:this.totNumShiftScen)';
                linearMask(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,3) = (1:this.totNumRangeScen)';
                
                [scenarios,ia] = unique(scenarios,'rows','stable');
                linearMask = linearMask(ia,:);
            
            else
                error('Not implemented yet!');
            end
                

            this.totNumScen = size(scenarios,1);


            this.relRangeShift = scenarios(:,5);
            this.absRangeShift = scenarios(:,4);
            this.isoShift = scenarios(:,1:3);
            
            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

            this.scenMask = false(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);
            
            this.scenForProb = scenarios;
            this.linearMask = linearMask;
            
            maskIx = sub2ind(size(this.scenMask),linearMask(:,1),linearMask(:,2),linearMask(:,3));
            this.scenMask(maskIx) = true;

            %Get Scenario probability
            Sigma = diag([this.shiftSD,this.rangeAbsSD,this.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            this.scenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios/cs).^2, 2)) / prod(diag(cs));
            this.scenWeight = this.scenProb./sum(this.scenProb);           
        end
    
        %% set methods
        function set.includeNominalScenario(this,value)
            arguments
                this
                value (1,1) {mustBeNumericOrLogical}
            end

            this.includeNominalScenario = value;
            this.update();
        end
    end


end

