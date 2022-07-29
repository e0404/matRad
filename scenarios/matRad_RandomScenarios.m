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
            
            %Scenario Probability from pdf
            this.scenForProb = scenarios;
            this.scenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios/cs).^2, 2)) / prod(diag(cs));

            %Scenario weight
            this.scenWeight = ones(this.nSamples,1)./this.nSamples; %equal weights since they have been randomly sampled (not entirely true if the Nominal scenario was forced) 

            %set variables
            this.totNumShiftScen = this.nSamples;
            this.totNumRangeScen = this.nSamples;
            this.totNumScen = this.nSamples; %check because of CT scenarios
            %this.totNumCtScen = 
            %this.numOfShiftScen = [nSamples,nSamples,nSamples];
            %this.numOfRangeShiftScen = nSamples;
            
            %Individual shifts
            this.relRangeShift = scenarios(:,5);
            this.absRangeShift = scenarios(:,4);
            this.isoShift = scenarios(:,1:3);

            this.maxAbsRangeShift = max(this.absRangeShift);
            this.maxRelRangeShift = max(this.relRangeShift);

            %Mask for scenario selection
            this.scenMask = false(this.numOfCtScen,this.totNumShiftScen,this.totNumRangeScen);

            for sCt = 1:this.numOfCtScen
                this.scenMask(sCt,:,:) = diag(true(this.nSamples,1));
            end
            
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            totNumScen    = sum(this.scenMask(:));

            if totNumScen ~= this.totNumScen
                matRad_cfg.dispWarning('Check Implementation of Total Scenario computation - given %d but found %d!',this.totNumScen,totNumScen);
                this.totNumScen = totNumScen;
            end
        end
    end


end

