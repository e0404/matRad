classdef matRad_NominalScenario < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements a single nominal planning scenario
%
% constructor
%   matRad_NominalScenario()
%   matRad_NominalScenario(ct)
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

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            this.updateScenarios();
        end
        
        function scenarios = updateScenarios(this)
            %Scenario Probability from pdf - here it is one since only one
            %scenario exist
            %TODO: In the context of an uncertainty model, we should
            %consider assigning probability according to the model, and
            %just leaving the weight 1
            this.scenForProb = zeros(this.numOfCtScen,5);
            this.scenProb = this.phaseProbability;

            %Scenario weight 
            this.scenWeight = ones(this.numOfCtScen,1)./this.numOfCtScen;
            this.scenWeight = this.phaseProbability;

            %set variables
            this.totNumShiftScen = 1;
            this.totNumRangeScen = 1;
            this.totNumScen = this.numOfCtScen; 
            
            %Individual shifts
            this.relRangeShift = zeros(this.numOfCtScen,1);
            this.absRangeShift = zeros(this.numOfCtScen,1);
            this.isoShift = zeros(this.numOfCtScen,3);

            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.absRangeShift);

            %Mask for scenario selection
            this.scenMask = true(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);
            
            %generic code
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            totNumScen    = sum(this.scenMask(:));
            
            %Get Scenario probability
            Sigma = diag([this.shiftSD,this.rangeAbsSD,this.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            tmpScenProb = (2*pi)^(-d/2) * exp(-0.5*sum((this.scenForProb/cs).^2, 2)) / prod(diag(cs));
            
            %Multiply with 4D phase probability
            this.scenProb = this.phaseProbability .* tmpScenProb;
            
            %Get relative (normalized) weight of the scenario
            this.scenWeight = this.scenProb./sum(this.scenProb); 
            
            %Return variable
            scenarios = (1:this.numOfCtScen)';
            scenarios = [scenarios zeros(this.numOfCtScen,5)];
            this.scenForProb = scenarios;

            if totNumScen ~= this.totNumScen
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
        
    end
end

