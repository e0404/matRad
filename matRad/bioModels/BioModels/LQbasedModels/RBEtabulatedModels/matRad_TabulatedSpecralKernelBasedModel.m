classdef matRad_TabulatedSpecralKernelBasedModel < matRad_LQRBETabulatedModel
    
    properties (Constant)
        model = 'TAB';
    end

    properties
        weightBy = 'Fluence';
        fragmentsToInclude = {'H'};
    end

    methods
        function this = matRad_TabulatedSpecralKernelBasedModel()
            this@matRad_LQRBETabulatedModel();
      
            % Cannot require energies and Data because these get
            % interpolated then by PB engine for example. PB takes the
            % quantities required here and interpolates them as kernels,
            % but they could also be something different. Need to find
            % different solution

            % this.requiredQuantities not assigned here because it depends
            % on the specific user defined properties. Could add here just
            % a default
            this.availableRadiationModalities = {'protons', 'carbon', 'helium'};
        end


        function bixel = calcBiologicalQuantitiesForBixel(this,bixel,kernel)
        
            bixel = calcBiologicalQuantitiesForBixel@matRad_LQRBETabulatedModel(this,bixel);

            nFragments = numel(this.fragmentsToInclude);    
            bixel.spectraEnergies = cellfun(@(fragment) bixel.baseData.Spectra.(this.weightBy).(fragment).energies,this.fragmentsToInclude, 'UniformOutput',false);
            
            % Get the tissue classes within the bixel
            bixelTissueIndexes = unique(bixel.vTissueIndex)';
            
            % Interpolate the alpha/beta table for the specific fragment (including primaries)
            % and alphaX/betaX ratio (tissue class)
            for i=bixelTissueIndexes
                [alphaE(i,:), betaE(i,:)] = arrayfun(@(fragment) this.interpolateRBETableForBixel(bixel.spectraEnergies{fragment}, this.fragmentsToInclude{fragment}, i),[1:nFragments], 'UniformOutput',false);
            end


            % Get the spectra kernels to be used (one for each fragment).
            kernelName = arrayfun(@(fragment) this.requiredQuantities{fragment}(this.requiredQuantities{fragment} ~= '.'),[1:nFragments], 'UniformOutput', false);

            bixelSpectra = arrayfun(@(fragment) kernel.(kernelName{fragment}),[1:nFragments], 'UniformOutput', false);

           
            % Get normalization for each fragment
            spectraFragmentDenominator = arrayfun(@(fragment) sum(bixelSpectra{fragment},2),[1:nFragments], 'UniformOutput',false);

            % Get total normalization by summing over the fragments
            spectraBixelDenominator = sum([spectraFragmentDenominator{:}],2);

            % Create containers for alpha and beta for each fragment
            alphaWeightedFragmentSpectra = NaN*ones(numel(bixel.radDepths),nFragments);
            betaWeightedFragmentSpectra = NaN*ones(numel(bixel.radDepths),nFragments);
            
            for i=bixelTissueIndexes
                % Get the spectrum-weighted alpha/beta for each fragment (sums over the energy bins at each bixel radDepth)
                tmpAlphaWeightedFragmentSpectra = arrayfun(@(fragment) bixelSpectra{fragment}*alphaE{i,fragment}, [1:nFragments], 'UniformOutput',false);
                tmpBetaWeightedFragmentSpectra  = arrayfun(@(fragment) bixelSpectra{fragment}*betaE{i,fragment}, [1:nFragments], 'UniformOutput',false);
                
                % Put it in matrix form and sum over the fragments. Only
                % include values for the correct tissue class for each
                % radDepth
                matTmpAlphaWeightedFragmentSpectra = [tmpAlphaWeightedFragmentSpectra{:}];
                matTmpBetaWeightedFragmentSpectra  = [tmpBetaWeightedFragmentSpectra{:}];

                alphaWeightedSpectra(bixel.vTissueIndex == i,:) = sum(matTmpAlphaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
                betaWeightedSpectra(bixel.vTissueIndex == i,:)  = sum(matTmpBetaWeightedFragmentSpectra(bixel.vTissueIndex == i,:),2);
            end

            % Normalize the quantities
            bixel.alpha = alphaWeightedSpectra./spectraBixelDenominator;
            bixel.beta  = betaWeightedSpectra./spectraBixelDenominator;
        end


        function assignBioModelPropertiesFromEngine(this, engine)

            matRad_cfg = MatRad_Config.instance();

            % Call superclass funtion
            assignBioModelPropertiesFromEngine@matRad_LQRBETabulatedModel(this, engine);

            % Check fragments
            if isprop(engine, 'bioProperties') 
                
                if isfield(engine.bioProperties, 'fragmentsToInclude')
                    this.fragmentsToInclude = engine.bioProperties.fragmentsToInclude;
                else
                    matRad_cfg.dispWarning('No fragments included! Only using ions identical to primary, this might result in inaccurate prediciton!');
                    switch engine.machine.meta.radiationMode
                        case 'protons'
                            this.fragmentsToInclude = {'H'};
                        case 'carbon'
                            this.fragmentsToInclude = {'C'};
                        case 'helium'
                            this.fragmentsToInclude = {'He'};
                    end
                end

                % Check the weighting factors (spectra)
                if  isfield(engine.bioProperties, 'weightBy')
                    this.weightBy = engine.bioProperties.weightBy;
                else
                    this.weightBy = 'Fluence';
                end

            end

            % This field is checked in the base data, depends on the
            % specific spectra and fragments that are included by the user
            % This also becomes the name tag for the bixel field containing
            % this interpolated quantity
            this.requiredQuantities = cellfun(@(fragment) ['Spectra.', this.weightBy, '.', fragment, '.Data'], this.fragmentsToInclude, 'UniformOutput',false); %[requiredSpectraData;requriedEnergies];

            % Check table consistency
            this.checkTableConsistency();
        end


        function checkTableConsistency(this)

            matRad_cfg = MatRad_Config.instance();
            % For the time being this checks for the RBEtabel field
            % data.includedIonZ. TODO: change in the RBE table entry the
            % included Ion Z with some more structured particle information
            % nABratios = numel(this.RBEtable.data);


            %availableZs = arrayfun(@(abRatio) abRatio.includedIonZ, this.RBEtable.data, 'UniformOutput',false);
            availableZs = this.RBEtable.data(1).includedIons;
            requiredFragments = this.fragmentsToInclude;
            
            if any(~ismember(this.fragmentsToInclude, availableZs))
                excludedFragments = this.fragmentsToInclude(~ismember(this.fragmentsToInclude, availableZs));
                matRad_cfg.dispError('Included RBE table does not contain information %s,',excludedFragments);
            end
          
        end
    end
end