classdef (Abstract) matRad_LQRBETabulatedModel < matRad_LQBasedModel
% This is an Abstract class implementig a tabulated RBE model.
% The model can handle multiple tissue alphaX/betaX ratio specified by the 
% cst structure, as long as a compatible RBEtable is provided.
%
% Properties of the model that can be defined by the user in pln.propDoseCalc.bioProperties:
%   RBEtable:       name of the specific table to be included
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    properties
        RBEtableName;
        defaultRBETable;
        fragmentsToInclude;
    end

    properties (SetAccess = protected, GetAccess = public)
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
            %matRad_cfg = MatRad_Config.instance();
            
            this@matRad_LQBasedModel();
            
            this.assignDefaultProperties();
            
            % This just for testing
            %this.defaultRBETable = 'RBEtable_rapidLEM_Russo2011_longErange_LEMI30';
        
        end

        function [alphaE,betaE] = interpolateRBETableForBixel(this,interpEnergies, fragment, tissueClass)
            % This function interpolates the correct table for the input
            % fragment and tissue class
            
            fragmentRBEtable = this.selectDataTableForFragment(fragment,tissueClass);

            alphaE = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.alpha, interpEnergies);
            betaE  = matRad_interp1(fragmentRBEtable.energies, fragmentRBEtable.beta,  interpEnergies);

        end

        function vTissueIndex = getTissueInformation(this,~, cst, dij,vAlphaX, ~, VdoseGrid, VdoseGridScenIx)

            % This function assumes that info in cst{i,5} has Tissue class =
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
            
            % This function selects the specific fragment and tissue class
            % entry in the RBE table.
            % The RBE table should be a struct containing the following
            % fields:
            %   .meta   struct with meta information
            %   .data   struct array with one entry for each tissue class
            % In turn, the data struct should contain the subfields:
            %   includedIons    ions included for the specific table (i.e 'H', 'He', 'C', ...);                                      
            %   energies        array containing the energies corresponding
            %                   to the specified alpha/beta tables
            %   alpha           matrix (#energies,#fragments) specifieng the alpha
            %                   value for the included fragments and
            %                   energies
            %   beta            matrix (#energies,#fragments) specifieng the alpha
            %                   value for the included fragments and
            %                   energies

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

    methods

        function assignDefaultProperties(this)
            this.defaultRBETable =  'RBEtable_rapidLEMI_testTable';
            this.fragmentsToInclude = {'H'};
        end

        function updateRBEtable(this)
            
            this.RBEtable = this.loadRBEtable(this.RBEtableName);

            this.availableAlphaInTable = [this.RBEtable.data(:).alphaX]';
            this.availableBetaInTable  = [this.RBEtable.data(:).betaX]';
            this.availableFragmentsInTable = [this.RBEtable.data(1).includedIons];
        end

      
    end

    
    methods %(Setters)
        
        function set.RBEtableName(this, value)

            this.RBEtableName = value;
            this.updateRBEtable();

        end

        function set.fragmentsToInclude(this, value)

            if iscell(value)
                this.fragmentsToInclude = value;
                this.updatePropertyValues();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Fragments to include should be a cell array: %s', tostring(value));
            end

        end
    
    end

    methods (Static)
        function RBEtable = loadRBEtable(fileName)
            
            % This function loads the specified RBE table
            matRad_cfg = MatRad_Config.instance();

            searchPath = {fullfile(matRad_cfg.matRadSrcRoot,'bioModels','RBEtables'),...    % default matratd folder
                          fullfile(matRad_cfg.primaryUserFolder, 'RBEtables')};             % user defined RBE table

            try
                load(fullfile(searchPath{1}, [fileName, '.mat']), 'RBEtable');
            catch
                try
                    laod(fullfile(searchPath{2}, [fileName, '.mat']), 'RBEtable');
                catch
                    matRad_cfg.dispError('Cannot find RBEtable: %s', fileName);
                end
            end

        end

        function checkTableConsistency(RBEtable, fragments)

           
            matRad_cfg = MatRad_Config.instance();
            % For the time being this checks for the RBEtabel field
            % data.includedIonZ. TODO: change in the RBE table entry the
            % included Ion Z with some more structured particle information
            % nABratios = numel(this.RBEtable.data);


            %availableZs = arrayfun(@(abRatio) abRatio.includedIonZ, this.RBEtable.data, 'UniformOutput',false);
            availableZs = RBEtable.data(1).includedIons;
            %requiredFragments = this.fragmentsToInclude;
            
            if any(~ismember(fragments, availableZs))
                excludedFragments = fragments(~ismember(fragments, availableZs));
                matRad_cfg.dispError('Included RBE table does not contain information %s,',excludedFragments);
            end
          
        end


        function checkRBEtableStructure(RBEtable)
            % Additional function to check structure of the Table

            matRad_cfg = MatRad_Config.instance();

            % Check structure of the table
            if ~isstruct(RBEtable)
                matRad_cfg.dispError('Provided Table is not a struct!');
            end

            if ~isequal(fieldnames(RBEtable), {'meta', 'data'}')

                matRad_cfg.dispError('Provided Table does not contain meta and data fields');

            end

            numOfTissues = numel(RBEtable.data);
            tissueAlphaBetaRatios = [RBEtable.data.alphaX]./[RBEtable.data.betaX];

            includedIons = RBEtable.data(1).includedIons;

            matRad_cfg.dispInfo('%2u tissues found in RBEtable\n', numOfTissues);
            matRad_cfg.dispInfo('with alpha/beta reatios of:\n');

            for tIdx = tissueAlphaBetaRatios
                matRad_cfg.dispInfo('\t %2.3f:\n', tIdx);
            end

            matRad_cfg.dispInfo('Data available for ions: ');
            for fIdx = includedIons
                matRad_cfg.dispInfo('%s ', fIdx{1});
            end

            matRad_cfg.dispInfo('\n');

        end
    end
end