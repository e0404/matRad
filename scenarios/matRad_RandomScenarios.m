classdef matRad_RandomScenarios < matRad_ScenarioModel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        includeNominalScenario = false;
        name = 'rndScen';
    end

    properties (Access = private)
        nSamplesDefault = 10;
    end
    
    methods
        function this = matRad_RandomScenarios(ct)           
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_ScenarioModel(superclassArgs{:});
        end
        
        function scenarios = updateScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            nSamples = this.nSamplesDefault;
            
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
            scenarios = randn(nSamples,d) * cs;

            if this.includeNominalScenario
                %We include the nominal scenario by just replacing the  
                %first one to keep the number of scenarios the same 
                scenarios(1,:) = 0; 
            end
            
            %Scenario Probability from pdf
            this.scenForProb = scenarios;
            this.scenProb = (2*pi)^(-d/2) * exp(-0.5*sum((scenarios/cs).^2, 2)) / prod(diag(cs));

            %Scenario weight
            this.scenWeight = ones(nSamples,1)./nSamples; %equal weights since they have been randomly sampled (not entirely true if the Nominal scenario was forced) 

            %set variables
            this.totNumShiftScen = nSamples;
            this.totNumRangeScen = nSamples;
            this.totNumScen = nSamples; %check because of CT scenarios
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
                this.scenMask(sCt,:,:) = diag(true(nSamples,1));
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

