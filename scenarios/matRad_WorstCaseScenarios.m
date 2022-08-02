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
            matRad_cfg = MatRad_Config.instance();
            
            %Create worst case setup shifts
            wcSetupShifts = this.wcSigma * this.shiftSD;            
            switch this.combinations
                case 'none'
                    %Create independent shifts
                    wcSetupShifts = wcSetupShifts.*[eye(3); -eye(3)];
                case 'shift'
                    %error('Not implemented!');
                    if this.includeNominalScenario
                        wcSetupShifts = [wcSetupShifts; zeros(1,numel(wcSetupShifts))];
                    end
                    wcSetupShifts = [wcSetupShifts; -wcSetupShifts(1,:)];
                    
                    [X,Y,Z] = meshgrid(wcSetupShifts(:,1),wcSetupShifts(:,2),wcSetupShifts(:,3));
                    wcSetupShifts = [X(:), Y(:), Z(:)];

                case 'all'
                    error('Not implemented!');
            end
            this.totNumShiftScen = size(wcSetupShifts,1);
                                
            %Create worst case range shifts
            wcRangeShifts = this.wcSigma * [this.rangeAbsSD this.rangeRelSD];
            wcRangeShifts(2,:) = -wcRangeShifts;
            if ~this.combineRange
                wcRangeShifts(end+1,:) = this.wcSigma * [this.rangeAbsSD -this.rangeRelSD];
                wcRangeShifts(end+1,:) = this.wcSigma * [-this.rangeAbsSD this.rangeRelSD];
            end
            
            %consider percent
            wcRangeShifts(:,2) = wcRangeShifts(:,2)./100;                      
            
            %Aggregate scenarios
            this.totNumRangeScen = size(wcRangeShifts,1);

            if any(strcmp(this.combinations,{'none','shift'}))
                scenarios = zeros(this.totNumShiftScen + this.totNumRangeScen,5);
                scenarios(1:this.totNumShiftScen,1:3) = wcSetupShifts;
                scenarios(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,4:5) = wcRangeShifts;
                
                %create the linear mask of scenarios
                linearMask = ones(size(scenarios,1),3);
                if strcmp(this.combinations,'none')
                    linearMask(1:this.totNumShiftScen,2) = (1:this.totNumShiftScen)' + this.includeNominalScenario;
                    linearMask(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,3) = (1:this.totNumRangeScen)' + this.includeNominalScenario;
                elseif strcmp(this.combinations,'shift')
                    linearMask(1:this.totNumShiftScen,2) = (1:this.totNumShiftScen)' + this.includeNominalScenario;
                    linearMask(this.totNumShiftScen+1:this.totNumShiftScen+this.totNumRangeScen,3) = (1:this.totNumRangeScen)' + this.includeNominalScenario;
                else
                    matRad_cfg.dispError('This sanity check should never be reached!');
                end
            
            else
                error('Not implemented yet!');
            end
                

            if this.includeNominalScenario
                %We include the nominal scenario by just replacing the  
                %first one to keep the number of scenarios the same 
                scenarios = [0 0 0 0 0; scenarios];
                linearMask = [1 1 1; linearMask];
                this.totNumRangeScen = this.totNumRangeScen + 1;
                this.totNumShiftScen = this.totNumShiftScen + 1;
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
    end


end
