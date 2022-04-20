classdef matRad_ScenarioModel < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Uncertainty model
        rangeRelSD  = 3.5;                % given in %
        rangeAbsSD  = 1;                  % given in [mm]
        shiftSD     = [2.25 2.25 2.25];   % given in [mm]
        wcSigma     = 1;                  % Multiplier to compute the worst case / maximum shifts    

        %scenProbMode = 'probBins'; %'pointwise','equalProb'
        %probDist = 'normal';
    end
   
    properties (SetAccess = protected)
        numOfCtScen;           % total number of CT scenarios


        % these parameters will be filled according to the choosen scenario type
        isoShift;
        relRangeShift;
        absRangeShift;

        maxAbsRangeShift;
        maxRelRangeShift;
        
        totNumShiftScen;        % total number of shift scenarios in x,y and z direction
        totNumRangeScen;        % total number of range and absolute range scenarios
        totNumScen;             % total number of samples 
        
        scenForProb;            % matrix for probability calculation - each row denotes one scenario
        scenProb;               % probability of each scenario stored in a vector (according to uncertainty model)
        scenWeight;             % weight of scenario relative to the underlying uncertainty model (depends on how scenarios are chosen / sampled)
        scenMask;
        linearMask;
    end
    
    methods
        function this = matRad_ScenarioModel(ct)
            if nargin == 0 || isempty(ct)
                this.numOfCtScen = 1;
            else
                this.numOfCtScen = ct.numOfCtScen;
            end

            this.updateScenarios();
        end     

        function scenarios = updateScenarios(this)            
        end

        function newInstance = extractSingleScenario(this,scenIdx)

            newInstance = matRad_NominalScenario();
           


            %%%%%%%%%%%%%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%
            newInstance = this;
            newInstance.TYPE = 'nomScen';
            newInstance.lockInit = true;
            newInstance.lockUpdate = true;
                        
            newInstance.scenForProb         = this.scenForProb(i,:);
            newInstance.relRangeShift       = this.scenForProb(i,5);
            newInstance.absRangeShift       = this.scenForProb(i,4);
            newInstance.isoShift            = this.scenForProb(i,1:3);
            newInstance.totNumShiftScen     = 1;
            newInstance.totNumRangeScen     = 1;
            newInstance.numOfCtScen         = ctScen;
            newInstance.scenMask            = 1;
            newInstance.linearMask          = 1;
            newInstance.scenProb            = 1;
            
            newInstance.lockInit = false;
            newInstance.lockUpdate = false;
        end

        %function calcScenProb(this,samplePos)
    end

    methods (Static)
        
    end
end

