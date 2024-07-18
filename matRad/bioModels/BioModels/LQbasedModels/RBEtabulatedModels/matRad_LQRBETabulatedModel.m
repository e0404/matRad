classdef (Abstract) matRad_LQRBETabulatedModel < matRad_LQBasedModel

    properties
        RBEtableName;
        defaultRBETable;
        RBEtable;
    end

    properties (Hidden)
        tissueAlphaX;               % array containing the alphaX values given in cst{i,5}. Dimension is (1,max(TissueClasses))
        tissueBetaX;                % array containing the betaX  values given in cst{i,5}. Dimension is (1,max(TissueClasses))
        availableAlphaInTable;
        availableBetaInTable;
        availableFragmentsInTable;
    end

    methods
        function this = matRad_LQRBETabulatedModel()
            matRad_cfg = MatRad_Config.instance();
            
            this@matRad_LQBasedModel();
            % This just for testing
            this.defaultRBETable = 'RBEtable_rapidLEM_Russo2011_longErange_LEMI30';
        
        end

        function [alphaE,betaE] = interpolateRBETableForBixel(this,interpEnergies, fragment, tissueClass)
            
            fragmentRBEtable = this.selectDataTableForFragment(fragment,tissueClass);

            alphaE = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.alpha, interpEnergies);
            betaE  = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.beta,  interpEnergies);

        end

        function assignBioModelPropertiesFromEngine(this, engine)

            matRad_cfg = MatRad_Config.instance();

            % Check RBE table
            if isprop(engine, 'bioProperties')
                if isfield(engine.bioProperties, 'RBEtable')
                
                    this.RBEtableName = engine.bioProperties.RBEtable;
                else
                    this.RBEtableName = this.defaultRBETable;
                    matRad_cfg.dispWarning('No RBE table specified, using default table: %s', this.defaultRBETable)
                end

            end
           
            this.RBEtable = this.loadRBEtable(this.RBEtableName);

            this.availableAlphaInTable = [this.RBEtable.data(:).alphaX]';
            this.availableBetaInTable  = [this.RBEtable.data(:).betaX]';
            this.availableFragmentsInTable = [this.RBEtable.data(1).includedIons];

        end

        function vTissueIndex = getTissueInformation(this,~, cst, dij,vAlphaX, ~, VdoseGrid, VdoseGridScenIx)

            % This assumes that info in cst{i,5} has Tissue class =
            % integer, alphaX,betaX and combination of alphaX,betaX for
            % given tissue class is consistent
            matRad_cfg = MatRad_Config.instance();

            numOfCtScen = numel(vAlphaX);

            cstDownsampled = matRad_setOverlapPriorities(cst);

            % resizing cst to dose cube resolution
            cstDownsampled = matRad_resizeCstToGrid(cstDownsampled,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
            
            tmpScenVdoseGrid = cell(numOfCtScen,1);

            for s = 1:numOfCtScen            
                tmpScenVdoseGrid{s} = VdoseGrid(VdoseGridScenIx{s});
                vTissueIndex{s}     = zeros(size(tmpScenVdoseGrid{s},1),1);
            end

            allTissueClasses = cellfun(@(cstStruct) cstStruct.TissueClass, cstDownsampled(:,5), 'UniformOutput',false);

            uniqueTissueClasses = unique([allTissueClasses{:}]);

            this.tissueAlphaX = NaN*ones(numel(uniqueTissueClasses),1);
            this.tissueBetaX  = NaN*ones(numel(uniqueTissueClasses),1);

            for i=uniqueTissueClasses
                %These should correspond to the tissue classes sampled by
                %the bixel
                this.tissueAlphaX(i) =  cstDownsampled{find([allTissueClasses{:}]==i,1,'first'),5}.alphaX;
                this.tissueBetaX(i)  =  cstDownsampled{find([allTissueClasses{:}]==i,1,'first'),5}.betaX;
            end

            for i = 1:size(cstDownsampled,1)
                % check if cst is compatiable
                if ~isempty(cstDownsampled{i,5}) && isfield(cstDownsampled{i,5},'alphaX') && isfield(cstDownsampled{i,5},'betaX')
 
                    %Set mapping between alphaX,betaX and tissueIndex
                    
                    IdxTissue = cstDownsampled{i,5}.TissueClass;
 
                    % check consitency of biological baseData and cst settings
                    if ~isempty(IdxTissue)
                        for s = 1:numOfCtScen
                            tmpScenVdoseGrid = VdoseGrid(VdoseGridScenIx{s});
                            isInVdoseGrid = ismember(tmpScenVdoseGrid,cstDownsampled{i,4}{s});
                            vTissueIndex{s}(isInVdoseGrid) = IdxTissue;
                        end
                    end
 
                else
                    for s = 1:numOfCtScen
                        vTissueIndex{s}(:) = 1;
                    end
                    matRad_cfg.dispWarning('\tTissue type of %s was set to 1\n',cstDownsampled{i,2});
                end
            end

            
            matRad_cfg.dispInfo('done.\n');
            
        end

        function fragmentRBEtable = selectDataTableForFragment(this,fragment,tissueClass)
            
            matRad_cfg = MatRad_Config.instance();

            selectedAlphaX = this.tissueAlphaX(tissueClass);
            selectedBetaX  = this.tissueBetaX(tissueClass);

            ratioIndexInTable    = find(ismember([this.availableAlphaInTable, this.availableBetaInTable], [selectedAlphaX, selectedBetaX], 'rows'));
            fragmentIndexInTable = find(strcmp(fragment,this.availableFragmentsInTable));
            
            if ~isempty(fragmentIndexInTable)
                if ~isempty(ratioIndexInTable)
                    fragmentRBEtable.alpha    = this.RBEtable.data(ratioIndexInTable).alpha(:,fragmentIndexInTable);
                    fragmentRBEtable.beta     = this.RBEtable.data(ratioIndexInTable).beta(:,fragmentIndexInTable);
                    fragmentRBEtable.energies = this.RBEtable.data(ratioIndexInTable).energies;
                else
                    matRad_cfg.dispError('AlphaX/BetaX ratio = %f/%f not available in RBEtable: %s',selectedAlphaX,selectedBetaX,this.RBEtableName);
                end
            else
                matRad_cfg.dispError('fragment %s not available in table: %s', fragment, this.RBEtableName);
            end
        end

    end

    methods (Static)
        function RBEtable = loadRBEtable(fileName)
            %get available tables in folder
                % RBE table should have meta and data, data.(ParticleName),
                % and each particle has alpha/beta/energy/mass/charge
    
            % Should load the whole RBE table, this should have multiple data

            % entries for different alphaX/betaX
            matRad_cfg = MatRad_Config.instance();


            load(fullfile(matRad_cfg.matRadSrcRoot,'bioModels','RBEtables', [fileName, '.mat']), 'RBEtable');

        end


    end
end