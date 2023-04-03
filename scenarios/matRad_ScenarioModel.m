classdef matRad_ScenarioModel < handle
%  matRad_ScenarioModel
%  This is an abstract interface class to define Scenario Models for use in
%  robust treatment planning and uncertainty analysis.
%  Subclasses should at least implement the update() function to generate
%  their own scenarios.
%
% constructor (Abstract)
%   matRad_ScenarioModel()
%   matRad_ScenarioModel(ct)
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
    
    properties (AbortSet = true) %We use AbortSet = true here to avoid updates when 
        %Uncertainty model
        rangeRelSD  = 3.5;                % given in %
        rangeAbsSD  = 1;                  % given in [mm]
        shiftSD     = [2.25 2.25 2.25];   % given in [mm]
        wcSigma     = 1;                  % Multiplier to compute the worst case / maximum shifts    
    end

    properties (Abstract,SetAccess=protected)
        name
    end

    properties (Dependent)
        wcFactor;
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
            
            %TODO: We could do this here automatically in the constructur, but
            %Octave 5 has a bug here and throws an error
            %this.updateScenarios();
        end

        function listAllScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Listing all scenarios...\n');
            matRad_cfg.dispInfo('\t#\txShift\tyShift\tzShift\tabsRng\trelRng\tprob.\n');
            for s = 1:size(this.scenForProb,1)
                str = num2str(this.scenForProb(s,:),'\t%.3f');
                matRad_cfg.dispInfo('\t%d\t%s\t%.3f\n',s,str,this.scenProb(s));
            end
        end

        %% SETTERS & UPDATE
        function set.rangeRelSD(this,rangeRelSD)
            valid = isnumeric(rangeRelSD) && isscalar(rangeRelSD) && rangeRelSD >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeRelSD! Needs to be a real positive scalar!');
            end
            this.rangeRelSD = rangeRelSD;
            this.updateScenarios();
        end

        function set.rangeAbsSD(this,rangeAbsSD)
            valid = isnumeric(rangeAbsSD) && isscalar(rangeAbsSD) && rangeAbsSD >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeAbsSD! Needs to be a real positive scalar!');
            end
            this.rangeAbsSD = rangeAbsSD;
            this.updateScenarios();
        end

        function set.shiftSD(this,shiftSD)
            valid = isnumeric(shiftSD) && isrow(shiftSD) && numel(shiftSD) == 3;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shiftSD! Needs to be 3-element numeric row vector!');
            end
            this.shiftSD = shiftSD;
            this.updateScenarios();
        end

        function set.wcSigma(this,wcSigma)
            valid = isnumeric(wcSigma) && isscalar(wcSigma) && wcSigma >= 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for wcSigma! Needs to be a real positive scalar!');
            end
            this.wcSigma = wcSigma;
            this.updateScenarios();
        end

           

        function scenarios = updateScenarios(this)            
            %This function will always update the scenarios given the
            %current property settings

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This abstract function needs to be implemented!');
        end

        function newInstance = extractSingleScenario(this,scenIdx)
            newInstance = matRad_NominalScenario();
                        
            newInstance.scenForProb         = this.scenForProb(scenIdx,:);
            newInstance.relRangeShift       = this.scenForProb(scenIdx,5);
            newInstance.absRangeShift       = this.scenForProb(scenIdx,4);
            newInstance.isoShift            = this.scenForProb(scenIdx,1:3);
            newInstance.scenProb            = this.scenProb(scenIdx);
            newInstance.scenWeight          = this.scenWeight(scenIdx);
            newInstance.numOfCtScen         = this.numOfCtScen;
        end
        
        %% Deprecated functions / properties
        function newInstance = extractSingleNomScen(this,~,scenIdx)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The function extractSingleNomScen of the scenario class will soon be deprecated! Use extractSingleScenario instead!');
            newInstance = this.extractSingleScenario(scenIdx);
        end

        function t = TYPE(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property TYPE of the scenario class will soon be deprecated!');
            t = this.name;
        end

        function value = get.wcFactor(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            value = this.wcSigma;
        end

        function set.wcFactor(this,value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            this.wcSigma = value;
        end


    end

    methods (Static)
        %{
        %TODO: implement automatic collection of available scenario classes
 
        function metaScenarioModels = getAvailableModels()
            matRad_cfg = MatRad_Config.instance();
            
            %Use the root folder and the scenarios folder only
            folders = {matRad_cfg.matRadRoot,mfilename("fullpath")};

            %
        end
        %}

        function types = AvailableScenCreationTYPE()
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The function/property AvailableScenarioCreationTYPE of the scenario class will soon be deprecated!');
            %Hardcoded for compatability with matRad_multScen
            types = {'nomScen','wcScen','impScen','rndScen'};
        end
    end
end

