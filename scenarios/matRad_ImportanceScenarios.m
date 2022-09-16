classdef matRad_ImportanceScenarios < matRad_ScenarioModel
%  matRad_ImportanceScenarios
%  Implements gridded importance scenarios, i.e., weighted according to a
%  probability distribution. It is not advised to create "combined" grids
%  with a large number of grid-points as the curse of dimensionality will
%  quickly break memory requirements when putting this in to dose influence
%  matrix computation.
%  
%
% constructor
%   matRad_ImportanceScenarios()
%   matRad_ImportanceScenarios(ct)
%
% input
%   ct:                 ct cube
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (AbortSet = true)
        includeNominalScenario = true;        
        combinations = 'none'; %Can be 'none', 'shift', 'all' to ontrol creation of worst case combinations 
        combineRange = true; %Wether to treat absolute & relative range as one shift or as separate scenarios

        numOfSetupGridPoints = 9;
        numOfRangeGridPoints = 9;
    end

    properties (SetAccess=protected)
        name = 'impScen';
    end

    properties (Constant)
        validCombinationTypes = {'all','none','shift'};
    end
    
    methods
        function this = matRad_ImportanceScenarios(ct)           
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            this.updateScenarios();
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
        
        function scenarios = updateScenarios(this)

            [scenarios,linearMask,this.totNumShiftScen,this.totNumRangeScen] = this.createGriddedScenarios(this,this.includeNominalScenario,this.combinations,this.combineRange,this.numOfSetupGridPoints,this.numOfRangeGridPoints);
                
            %Finalize meta information
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
    end

    methods (Static)
        function [scenarios,linearMask,totNumShiftScen,totNumRangeScen] = createGriddedScenarios(scenarioModel,includeNominalScenario,combinations,combineRange,numOfSetupGridPoints,numOfRangeGridPoints)
            matRad_cfg = MatRad_Config.instance();

            %Get the maximum, i.e., worst case shifts
            wcSetupShifts = scenarioModel.wcSigma * scenarioModel.shiftSD;
            
            %% Create gridded setup shifts
            %Create grid vectors for setup shifts
            setupShiftGrid = zeros(numOfSetupGridPoints,numel(wcSetupShifts));            
            if mod(numOfSetupGridPoints,2) == 0 && includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Setup Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end

            for i = 1:numel(wcSetupShifts)
                setupShiftGrid(:,i) = linspace(-wcSetupShifts(i),wcSetupShifts(i),numOfSetupGridPoints);
                if includeNominalScenario 
                      
                    [~,ix] = min(abs(setupShiftGrid(:,i)));
                    setupShiftGrid(ix,i) = 0;    
                end  
            end
            
            %Now create vector of all shifts for different combinatorial
            %settings
            switch combinations
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
            
            if ~isempty(shiftNomScenIx) || includeNominalScenario
                if ~isempty(shiftNomScenIx)
                    griddedSetupShifts(shiftNomScenIx,:) = [];
                end
                griddedSetupShifts = [0 0 0; griddedSetupShifts];
            end
                        
            totNumShiftScen = size(griddedSetupShifts,1);
                                
            %% Create gridded range shifts
            %Obtain worst case range shifts
            wcRangeShifts = scenarioModel.wcSigma * [scenarioModel.rangeAbsSD scenarioModel.rangeRelSD./100];        
            
            rangeShiftGrid = zeros(numOfRangeGridPoints,numel(wcRangeShifts));            
            if mod(numOfRangeGridPoints,2) == 0 && includeNominalScenario
                matRad_cfg.dispWarning('Obtaining Range Shifts: Including the nominal scenario with even number of grid points creates asymmetrical shifts!');
            end

            for i = 1:numel(wcRangeShifts)
                rangeShiftGrid(:,i) = linspace(-wcRangeShifts(i),wcRangeShifts(i),numOfRangeGridPoints);
                if includeNominalScenario 
                      
                    [~,ix] = min(abs(rangeShiftGrid(:,i)));
                    rangeShiftGrid(ix,i) = 0;    
                end  
            end

            if combineRange
                griddedRangeShifts = rangeShiftGrid;
            else                
                [rngAbs,rngRel] = meshgrid(rangeShiftGrid(:,1),rangeShiftGrid(:,2));
                griddedRangeShifts = [rngAbs(:),rngRel(:)];
            end

            %Remove duplicate scenarios and update number of shifts
            griddedRangeShifts = matRad_ImportanceScenarios.uniqueStableRowsCompat(griddedRangeShifts); 

            rangeNomScenIx = find(all(griddedRangeShifts == zeros(1,2),2));            
            
            if ~isempty(rangeNomScenIx) || includeNominalScenario
                if ~isempty(rangeNomScenIx)
                    griddedRangeShifts(rangeNomScenIx,:) = [];
                end
                griddedRangeShifts = [0 0; griddedRangeShifts];
            end

            totNumRangeScen = size(griddedRangeShifts,1);
            
            
            %Aggregate scenarios
            switch combinations
                case {'none','shift'}
                    scenarios = zeros(totNumShiftScen + totNumRangeScen,5);
                    scenarios(1:totNumShiftScen,1:3) = griddedSetupShifts;
                    scenarios(totNumShiftScen+1:totNumShiftScen+totNumRangeScen,4:5) = griddedRangeShifts;



                    %create the linear mask of scenarios
                    linearMask = ones(size(scenarios,1),3);
                    linearMask(1:totNumShiftScen,2) = (1:totNumShiftScen)';
                    linearMask(totNumShiftScen+1:totNumShiftScen+totNumRangeScen,3) = (1:totNumRangeScen)';

                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMask = linearMask(ia,:);

                case 'all'
                    %Prepare scenario matrix by replicating shifts
                    %with the number of range scenarios
                    scenarios = repmat(griddedSetupShifts,totNumRangeScen,1);
                    scenarios = [scenarios zeros(size(scenarios,1),2)];
                    
                    %create the linear mask of scenarios
                    linearMask = ones(size(scenarios,1),3);
                    for r = 1:totNumRangeScen                    
                        offset = (r-1)*totNumShiftScen;
                        ixR = (offset + 1):(offset + totNumShiftScen);
                        scenarios(ixR,4:5) = repmat(griddedRangeShifts(r,:),totNumShiftScen,1);

                        %Set linear mask
                        linearMask(ixR,2) = (1:totNumShiftScen)';
                        linearMask(ixR,3) = r;
                    end


                    %create the linear mask of scenarios
                    [scenarios,ia] = matRad_ImportanceScenarios.uniqueStableRowsCompat(scenarios);
                    linearMask = linearMask(ia,:);
                otherwise
                    matRad_cfg.dispError('Invalid value for combinations! This sanity check should never be reached!');
            end
        end
    
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

