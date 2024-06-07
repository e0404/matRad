classdef matRad_RandomScenarios < matRad_ScenarioModel
%  matRad_RandomScenarios
%  Implements randomly sampled scenarios
%
% constructor
%   matRad_RandomScenarios()
%   matRad_RandomScenarios(ct)
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
        includeNominalScenario = false; %Forces inclusion of the nominal scenario        
        nSamples = 10;                  %Standard number of random samples
    end

    %Deprecated Properties that were used
    properties (Dependent)
        numOfShiftScen;
        numOfRangeShiftScen;
    end

    properties (SetAccess=protected)
        name = 'rndScen';
    end
    
    methods
        function this = matRad_RandomScenarios(ct)           
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

        %% Setters & Update
        function set.nSamples(this,nSamples)
            valid = isnumeric(nSamples) && isscalar(nSamples) && mod(nSamples,1) == 0 && nSamples > 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for nSamples! Needs to be a positive integer!');
            end
            this.nSamples = nSamples;
            this.updateScenarios();
        end
        
        function set.numOfShiftScen(this,numOfShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            
            %That's for downwards compatibility
            if ~isscalar(numOfShiftScen)
                numOfShiftScen = unique(numOfShiftScen);
            end

            this.nSamples = numOfShiftScen;
        end

        function  value = get.numOfShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfShiftScen of the scenario class will soon be deprecated! Use nSamples instead');            
            value = this.nSamples;
        end

        function set.numOfRangeShiftScen(this,numOfRangeShiftScen)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');           
            this.nSamples = numOfRangeShiftScen;
        end

        function  value = get.numOfRangeShiftScen(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property numOfRangeShiftScen of the scenario class will soon be deprecated! Use nSamples instead');
            value = this.nSamples;
        end
        
        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();

            this.numOfCtScen = size(this.ctScenProb,1);
            
            %Multivariate Normal Sampling
            Sigma = diag([this.shiftSD,this.rangeAbsSD,this.rangeRelSD./100].^2);
            d = size(Sigma,1);
            [cs,p] = chol(Sigma);
            
            %The lower part is there but commented for completeness, we do
            %not need it since we know we that sigma is PSD
            %if p ~= 0 %for the positive semi-definite case, we don't check for negative definite
            %      %      [V,L] = eig(Sigma);
            %      cs = sqrt(L) * V';
            %end
            %transform normal samples, mean is always zero
            scenarios = randn(this.nSamples,d) * cs;

            if this.includeNominalScenario
                %We include the nominal scenario by just replacing the  
                %first one to keep the number of scenarios the same 
                scenarios(1,:) = 0; 
            end

            %Handle 4D phases
            phases = repmat(this.ctScenProb(:,1)',size(scenarios,1),1);
            phases = phases(:);
            scenarios = horzcat(phases, repmat(scenarios,[this.numOfCtScen 1]));
            this.ctScenIx = phases;
            
            %Scenario Probability from pdf
            this.scenForProb = scenarios;
            
            %Create phase and setup/range probabilities
            tmpScenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios(:,2:end)/cs).^2, 2)) / prod(diag(cs));
            
            %Now we combine with the 4D ct phase probabilities (multiply)
            tmpPhaseProb = arrayfun(@(phase) this.ctScenProb(find(this.ctScenProb(:,1) == phase),2),phases);
            
            %Finalize Probability
            this.scenProb = tmpPhaseProb .* tmpScenProb;
                
            %Scenario weight
            this.scenWeight = [];
            for sCt = 1:this.numOfCtScen
                %equal weights within a phase since they have been randomly sampled 
                %(not entirely true if the Nominal scenario was forced) 
                this.scenWeight = [this.scenWeight; this.ctScenProb(sCt,2) * ones(this.nSamples,1)./this.nSamples];
            end

            %set variables
            this.totNumShiftScen = this.nSamples;
            this.totNumRangeScen = this.nSamples;
            this.totNumScen = this.nSamples * this.numOfCtScen; %check because of CT scenarios
            %this.totNumCtScen = 
            %this.numOfShiftScen = [nSamples,nSamples,nSamples];
            %this.numOfRangeShiftScen = nSamples;
            
            %Individual shifts
            this.relRangeShift = scenarios(:,6);
            this.absRangeShift = scenarios(:,5);
            this.isoShift = scenarios(:,2:4);

            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

            %Mask for scenario selection
            this.scenMask = false(this.numOfAvailableCtScen,this.totNumShiftScen,this.totNumRangeScen);

            for sCt = 1:this.numOfCtScen
                scenIx = this.ctScenProb(sCt,1);
                this.scenMask(scenIx,:,:) = diag(true(this.nSamples,1));
            end
            

            tmpScenMask = permute(this.scenMask,[3 2 1]);
            [x{3}, x{2}, x{1}] = ind2sub(size(tmpScenMask),find(tmpScenMask));
            this.linearMask    = cell2mat(x);
            totNumScen    = sum(this.scenMask(:));

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
    end


end

