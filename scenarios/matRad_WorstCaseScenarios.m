classdef matRad_WorstCaseScenarios < matRad_ScenarioModel
%  matRad_WorstCaseScenarios
%  Implements single worst-case shifts per dimension.%  
%
% constructor
%   matRad_WorstCaseScenarios()
%   matRad_WorstCaseScenarios(ct)
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
    
    properties
        includeNominalScenario = true;
        combinations = 'none'; %Can be 'none', 'shift', 'all' to ontrol creation of worst case combinations 
        combineRange = true; %Wether to treat absolute & relative range as one shift or as separate scenarios
    end

    properties (SetAccess=protected)
        name = 'wcScen';
    end

    properties (Constant)
        validCombinationTypes = {'all','none','shift'};
    end
    
    methods
        function this = matRad_WorstCaseScenarios(ct)           
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

        function set.includeNominalScenario(this,includeNomScen)
            valid = isscalar(includeNomScen) && (isnumeric(includeNomScen) || islogical(includeNomScen));
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for includeNominalScenario! Needs to be a boolean / logical value!');
            end
            this.includeNominalScenario = includeNomScen;
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
            
            if this.includeNominalScenario
                nPoints = 3;
            else
                nPoints = 2;
            end
            
            %Use the static gridded shift function from
            %ImportanceScenarios. We set inclusion of nominal scenarios to
            %false and handle it automatically via the grid point number
            [scenarios,linearMask,this.totNumShiftScen,this.totNumRangeScen] = matRad_ImportanceScenarios.createGriddedScenarios(this,false,this.combinations,this.combineRange,nPoints,nPoints);
            
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
    end


end
