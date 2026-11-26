classdef (Abstract) matRad_ScenarioModel < handle
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
    % input:
    %   ct:                 ct cube
    %
    % Scenario indexing:
    %   A scenario is addressed in three ways.
    %
    %   - CT scenario id: index of the CT scenario in ``ct.cubeHU`` (and
    %     ``cst{:,4}``), as listed in the first column of ``ctScenProb``.
    %     If only a subset of the CT scenarios is used, the CT scenario id
    %     differs from the position of the CT scenario in ``ctScenProb``.
    %   - Scenario number: row in ``linearMask``, ``scenForProb``,
    %     ``scenProb``, ``scenWeight``, ``isoShift``, ``absRangeShift`` and
    %     ``relRangeShift`` (``1:totNumScen``).
    %   - Full scenario index: linear index into ``scenMask``, which is also
    %     the index of the scenario in the dose influence cell arrays
    %     (e.g. ``dij.physicalDose``). Its subscripts are the CT scenario id,
    %     the shift scenario and the range shift scenario stored in
    %     ``linearMask``.
    %
    %   ``sub2scenIx`` and ``scenNum2scenIx`` convert to the full scenario
    %   index, ``scenNum`` converts the full scenario index back to the
    %   scenario number.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022-2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (AbortSet = true) % We use AbortSet = true here to avoid updates when
        % Uncertainty model
        rangeRelSD  = 3.5                 % given in %
        rangeAbsSD  = 1                   % given in [mm]
        shiftSD     = [2.25 2.25 2.25]    % given in [mm]
        wcSigma     = 1                   % Multiplier to compute the worst case / maximum shifts

        ctScenProb  = [1 1]               % Ct Scenarios to be included in the model. Left column: Scenario Index. Right column: Scenario Probability
    end

    properties (Abstract, SetAccess = protected)
        name
        shortName
    end

    properties (Dependent)
        wcFactor
    end

    properties (SetAccess = protected)
        numOfCtScen             % total number of CT scenarios used
        numOfAvailableCtScen    % total number of CT scenarios existing in ct structure
        ctScenIx                % map of all ct scenario indices per scenario
        motionPeriod = Inf      % motion period of 4D CT, if 4D CT

        % these parameters will be filled according to the chosen scenario type
        isoShift
        relRangeShift
        absRangeShift

        maxAbsRangeShift
        maxRelRangeShift

        totNumShiftScen         % total number of shift scenarios in x,y and z direction
        totNumRangeScen         % total number of range and absolute range scenarios
        totNumScen              % total number of samples

        scenForProb             % matrix for probability calculation - each row denotes one scenario, whereas columns denotes the realization value
        scenProb                % probability of each scenario stored in a vector (according to uncertainty model)
        scenWeight              % weight of scenario relative to the underlying uncertainty model (depends on how scenarios are chosen / sampled)
        scenMask
        linearMask
    end

    methods

        function this = matRad_ScenarioModel(ct)
            if nargin == 0 || isempty(ct)
                this.numOfCtScen = 1;
                this.numOfAvailableCtScen = 1;
            else
                this.numOfCtScen = ct.numOfCtScen;
                this.numOfAvailableCtScen = ct.numOfCtScen;
                if isfield(ct, 'motionPeriod')
                    this.motionPeriod  = ct.motionPeriod;
                end
            end

            % Equal probability to be in each phase of the 4D ct
            this.ctScenProb = [(1:this.numOfCtScen)', ones(this.numOfCtScen, 1) ./ this.numOfCtScen];

            % TODO: We could do this here automatically in the constructor, but
            % Octave 5 has a bug here and throws an error
            % this.updateScenarios();
        end

        function listAllScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Listing all scenarios...\n');
            matRad_cfg.dispInfo('\t#\tctScen\txShift\tyShift\tzShift\tabsRng\trelRng\tprob.\n');
            for s = 1:size(this.scenForProb, 1)
                str = [num2str(this.scenForProb(s, 1), '%d\t'), sprintf('\t'), num2str(this.scenForProb(s, 2:end), '\t%.3f')];
                matRad_cfg.dispInfo('\t%d\t%s\t%.3f\n', s, str, this.scenProb(s));
            end
        end

        %% SETTERS & UPDATE
        function set.rangeRelSD(this, rangeRelSD)
            valid = isnumeric(rangeRelSD) && isscalar(rangeRelSD) && rangeRelSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeRelSD! Needs to be a real positive scalar!');
            end
            this.rangeRelSD = rangeRelSD;
            this.updateScenarios();
        end

        function set.rangeAbsSD(this, rangeAbsSD)
            valid = isnumeric(rangeAbsSD) && isscalar(rangeAbsSD) && rangeAbsSD >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for rangeAbsSD! Needs to be a real positive scalar!');
            end
            this.rangeAbsSD = rangeAbsSD;
            this.updateScenarios();
        end

        function set.shiftSD(this, shiftSD)
            valid = isnumeric(shiftSD) && isrow(shiftSD) && numel(shiftSD) == 3 && all(shiftSD > 0);
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shiftSD! Needs to be 3-element numeric row vector!');
            end
            this.shiftSD = shiftSD;
            this.updateScenarios();
        end

        function set.wcSigma(this, wcSigma)
            valid = isnumeric(wcSigma) && isscalar(wcSigma) && wcSigma >= 0;
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for wcSigma! Needs to be a real positive scalar!');
            end
            this.wcSigma = wcSigma;
            this.updateScenarios();
        end

        function set.ctScenProb(this, ctScenProb)
            valid = isnumeric(ctScenProb) && ...
                ismatrix(ctScenProb) && ...
                size(ctScenProb, 2) == 2 && ...
                all(round(ctScenProb(:, 1)) == ctScenProb(:, 1)) && ...
                all(ctScenProb(:) >= 0);

            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError(['Invalid value for used ctScenProb! ' ...
                                      'Needs to be a valid 2-column matrix with left column representing the scenario index ' ...
                                      'and right column representing the appropriate probabilities [0,1]!']);
            end
            this.ctScenProb = ctScenProb;
            this.updateScenarios();
        end

        function scenarios = updateScenarios(this)
            % This function will always update the scenarios given the
            % current property settings

            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispError('This abstract function needs to be implemented!');
        end

        function newInstance = extractSingleScenario(this, scenNum)
            newInstance = matRad_NominalScenario();

            % The CT scenario id is not necessarily the row in ctScenProb
            ctScenId = this.linearMask(scenNum, 1);
            ctScenProbIx = find(this.ctScenProb(:, 1) == ctScenId, 1, 'first');
            if isempty(ctScenProbIx)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Could not find CT scenario %d in ctScenProb!', ctScenId);
            end

            % First set properties that force an update
            newInstance.numOfCtScen         = 1;
            newInstance.ctScenProb          = this.ctScenProb(ctScenProbIx, :);

            % Now overwrite existing variables for correct probabilities
            % and error realizations
            newInstance.scenForProb         = this.scenForProb(scenNum, :);
            newInstance.relRangeShift       = this.scenForProb(scenNum, 6);
            newInstance.absRangeShift       = this.scenForProb(scenNum, 5);
            newInstance.isoShift            = this.scenForProb(scenNum, 2:4);
            newInstance.scenProb            = this.scenProb(scenNum);
            newInstance.scenWeight          = this.scenWeight(scenNum);
            newInstance.maxAbsRangeShift    = max(abs(this.absRangeShift(scenNum)));
            newInstance.maxRelRangeShift    = max(abs(this.relRangeShift(scenNum)));
            newInstance.scenMask            = false(this.numOfAvailableCtScen, 1, 1);
            newInstance.linearMask          = [newInstance.ctScenIx 1 1];

            newInstance.scenMask(newInstance.linearMask(:, 1), newInstance.linearMask(:, 2), newInstance.linearMask(:, 3)) = true;
            % newInstance.updateScenarios();
        end

        function scenIx = sub2scenIx(this, ctScen, shiftScen, rangeShiftScen, ctScenReference)
            % Returns the full scenario index (linear index into scenMask
            % and the dose influence cell arrays) from scenario subscripts
            %
            % call:
            %   scenIx = model.sub2scenIx(ctScen, shiftScen, rangeShiftScen)
            %   scenIx = model.sub2scenIx(ctScen, shiftScen, rangeShiftScen, ctScenReference)
            %
            % input:
            %   ctScen:          CT scenario, interpreted according to
            %                    ctScenReference
            %   shiftScen:       shift scenario index
            %   rangeShiftScen:  range shift scenario index
            %   ctScenReference: (optional) 'id' (default) if ctScen is the
            %                    CT scenario id as stored in linearMask(:,1),
            %                    or 'position' if ctScen is the position of
            %                    the CT scenario in ctScenProb
            %
            % output:
            %   scenIx:          full scenario index

            matRad_cfg = MatRad_Config.instance();

            if nargin < 5 || isempty(ctScenReference)
                ctScenReference = 'id';
            end

            if ~ischar(ctScenReference) || ~any(strcmpi(ctScenReference, {'id', 'position'}))
                matRad_cfg.dispError('ctScenReference must be ''id'' or ''position''!');
            end

            if strcmpi(ctScenReference, 'position')
                if any(ctScen(:) < 1 | ctScen(:) > size(this.ctScenProb, 1) | mod(ctScen(:), 1) ~= 0)
                    matRad_cfg.dispError('CT scenario position exceeds the number of CT scenarios in ctScenProb!');
                end
                ctScen = this.ctScenProb(ctScen, 1);
            end

            maskSize = this.getScenMaskSize();
            subs = {ctScen, shiftScen, rangeShiftScen};
            for d = 1:3
                if any(subs{d}(:) < 1 | subs{d}(:) > maskSize(d) | mod(subs{d}(:), 1) ~= 0)
                    matRad_cfg.dispError('Scenario subscript exceeds the dimensions of the scenario mask!');
                end
            end

            scenIx = sub2ind(maskSize, ctScen, shiftScen, rangeShiftScen);
        end

        function scenIx = scenNum2scenIx(this, scenNum)
            % Returns the full scenario index (linear index into scenMask
            % and the dose influence cell arrays) of scenario numbers, i.e.,
            % rows of linearMask
            %
            % call:
            %   scenIx = model.scenNum2scenIx()
            %   scenIx = model.scenNum2scenIx(scenNum)
            %
            % input:
            %   scenNum:    (optional) scenario numbers, all scenarios if
            %               omitted
            %
            % output:
            %   scenIx:     column vector of full scenario indices

            if nargin < 2
                scenNum = 1:size(this.linearMask, 1);
            end

            subs = this.linearMask(scenNum, :);
            scenIx = sub2ind(this.getScenMaskSize(), subs(:, 1), subs(:, 2), subs(:, 3));
        end

        function scenNum = scenNum(this, fullScenIx)
            % Returns the scenario number, i.e., the row in linearMask,
            % scenForProb, scenProb, ..., from full scenario indices into
            % scenMask. Indices of inactive scenarios are dropped.
            %
            % call:
            %   scenNum = model.scenNum(fullScenIx)
            %
            % input:
            %   fullScenIx: full scenario indices
            %
            % output:
            %   scenNum:    scenario numbers

            [isActive, scenNum] = ismember(fullScenIx, this.scenNum2scenIx());
            scenNum = scenNum(isActive);
        end

        %% Deprecated functions / properties
        function newInstance = extractSingleNomScen(this, ~, scenIdx)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning(['The function extractSingleNomScen of the scenario class will soon be deprecated! ' ...
                                               'Use extractSingleScenario instead!']);
            newInstance = this.extractSingleScenario(scenIdx);
        end

        function t = TYPE(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property TYPE of the scenario class will soon be deprecated!');
            t = this.shortName;
        end

        function value = get.wcFactor(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            value = this.wcSigma;
        end

        function set.wcFactor(this, value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('The property wcFactor of the scenario class will soon be deprecated!');
            this.wcSigma = value;
        end

    end

    methods (Access = protected)

        function maskSize = getScenMaskSize(this)
            % Size of scenMask as [ctScen shiftScen rangeShiftScen], also
            % if trailing singleton dimensions are dropped by MATLAB
            maskSize = [size(this.scenMask) 1 1];
            maskSize = maskSize(1:3);
        end

    end

    methods (Static)

        function classList = getAvailableModels()
            matRad_cfg = MatRad_Config.instance();

            % Use the root folder and the scenarios folder only
            folders = {fileparts(mfilename('fullpath'))};
            folders = [folders matRad_cfg.userfolders];

            persistent metaScenarioModels lastOptionalPaths

            % First we do a sanity check if persistently stored metaclasses are valid
            if ~matRad_cfg.isOctave && ~isempty(metaScenarioModels) && ~all(cellfun(@isvalid, metaScenarioModels))
                matRad_cfg.dispWarning('Found invalid ScenarioModels, updating model cache.');
                metaScenarioModels = [];
            end

            if isempty(metaScenarioModels) || (~isempty(lastOptionalPaths) && ~isequal(lastOptionalPaths, folders))
                lastOptionalPaths = folders;
                metaScenarioModels = matRad_findSubclasses(meta.class.fromName(mfilename('class')), 'folders', folders, 'includeSubfolders', true);
            end
            classList = matRad_identifyClassesByConstantProperties(metaScenarioModels, 'shortName', 'defaults', {'nomScen'});

            if isempty(classList)
                matRad_cfg.dispError('No models found in paths %s', strjoin(folders, '\n'));
            end
        end

        function model = create(modelMetadata, ct)
            if isa(modelMetadata, 'matRad_ScenarioModel')
                model = modelMetadata;
                return
            end

            matRad_cfg = MatRad_Config.instance();

            if ischar(modelMetadata) || isstring(modelMetadata)
                modelMetadata = struct('model', modelMetadata);
            end

            modelClassList = matRad_ScenarioModel.getAvailableModels();
            modelNames = {modelClassList.shortName};

            if ~isfield(modelMetadata, 'model')
                matRad_cfg.dispWarning('Scenario Model not found, creating nominal scenario instead!');
                modelMetadata.model = 'nomScen';
            elseif ~(ischar(modelMetadata.model) || isstring(modelMetadata.model))
                matRad_cfg.dispSelectionFallback('Scenario model selection is not a model name, creating nominal scenario instead!');
                modelMetadata.model = 'nomScen';
            elseif ~any(strcmp(modelNames, modelMetadata.model))
                matRad_cfg.dispSelectionFallback('Scenario model ''%s'' not found (available: %s), creating nominal scenario instead!', ...
                                                 char(modelMetadata.model), strjoin(modelNames, ', '));
                modelMetadata.model = 'nomScen';
            end

            usedModel = find(strcmp(modelNames, modelMetadata.model));

            if ~isscalar(usedModel)
                usedModel = usedModel(1);
            end

            modelClassInfo = modelClassList(usedModel);

            if nargin < 2
                model = modelClassInfo.handle();
            else
                model = modelClassInfo.handle(ct);
            end

            modelMetadata = rmfield(modelMetadata, 'model');

            % Now overwrite properties
            fields = fieldnames(modelMetadata);

            % iterate over all fieldnames and try to set the
            % corresponding properties inside the engine
            if matRad_cfg.isOctave
                c2sWarningState = warning('off', 'Octave:classdef-to-struct');
            end

            for i = 1:length(fields)
                try
                    field = fields{i};
                    if matRad_ispropCompat(model, field)
                        model.(field) = matRad_recursiveFieldAssignment(model.(field), modelMetadata.(field), true);
                    else
                        matRad_cfg.dispWarning('Not able to assign property ''%s'' from multScen struct to Scenario Model!', field);
                    end
                catch ME
                    % catch exceptions when the model has no properties,
                    % which are defined in the struct.
                    % When defining an engine with custom setter and getter
                    % methods, custom exceptions can be caught here. Be
                    % careful with Octave exceptions!
                    matRad_cfg = MatRad_Config.instance();
                    switch ME.identifier
                        case 'MATLAB:noPublicFieldForClass'
                            matRad_cfg.dispWarning('Not able to assign property from multScen struct to scenario model: %s', ME.message);
                        otherwise
                            matRad_cfg.dispWarning('Problem while setting up scenario Model from struct:%s %s', field, ME.message);
                    end
                end
            end

            if matRad_cfg.isOctave
                warning(c2sWarningState.state, 'Octave:classdef-to-struct');
            end
        end

    end
end
