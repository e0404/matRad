classdef (Abstract) matRad_BiologicalModel < handle
    %  matRad_BiologicalModel
    %  This is an abstract interface class to define Biological Models for use in
    %  dose calculation and plan optimization.
    %  Subclasses should at least implement the methods:
    %
    %  calcBiologicalQuantitiesForBixel()      to implement the specific model calculation algorithm
    %
    %  Additional methods that can be implemented are:
    %
    %   getTissueInformation()                  to collect meta information about the tissue defintion and paramters
    %   assignBioModelPropertiesFromEngine()    to translate user defined paramters directly to the biological model subclass
    %
    %
    % All subclasses should also declare the  properties:
    %
    %   'possibleRadiationModes'         to specify the radiation modalities to which the model validity is limited
    %   'requiredQuantities'                   to check the availability of information stored in the provided machine file
    %
    % constructor (Abstract)
    %   matRad_BiologicalModel()
    %
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
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Abstract, Constant)
        model;
        requiredQuantities;     % kernels in base data needed for the alpha/beta calculation
        possibleRadiationModes; % radiation modalitites compatible with the model
        defaultReportQuantity;  % default suggested quantity to use for display and planning    
    end

    properties (Hidden)
        quantityOpt;
        quantityVis;
    end

    methods %(Abstract)
        function bixel = calcBiologicalQuantitiesForBixel(this)
            % the actual calculation method wich calculates biological quantities for individual beamlets
            % Needs to be implemented in non abstract subclasses. 
            throw(MException('MATLAB:class:AbstractMember','Abstract function calcBiologicalQuantitiesForBixel of your BiologicalModel needs to be implemented!'));
        end
    end

    methods
        function set.quantityOpt(this, value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('Property quantityOpt is deprecated from bioModel. Please set it as a field in pln.propOpt');
        end

        function value = get.quantityOpt(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('Property quantityOpt is deprecated from bioModel.');
            value = [];
        end

        function set.quantityVis(this, value)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('Property quantityVis is deprecated from bioModel. Please set it as a field in pln.propOpt');
        end

        function value = get.quantityVis(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispDeprecationWarning('Property quantityVis is deprecated from bioModel.');
            value = [];
        end
    end

    methods %(Static)
        function [vTissueIndex] = getTissueInformation(this,~,~,~,vAlphaX,~,~,~) %(machine,cst,dij,vAlphaX,vBetaX,VdoseGrid, VdoseGridScenIdx)
            % This is the default, should be masked by the specific model
            % subclass if needed

            for s=1:numel(vAlphaX)
                vTissueIndex{s} = zeros(size(vAlphaX{s}));
            end
        end

        function [avail, msg] = isAvailable(this, radiationMode, providedQuantities)
            validRadMode = any(strcmp(this.possibleRadiationModes,radiationMode));
            
            msg1 = '';
            msg2 = '';
            if ~validRadMode
                msg1 = sprintf('Radiation mode invalid for model %s!',this.model);
            end

            % Quantity check
            validQuantities = true;
            if nargin > 2 && ~isempty(this.requiredQuantities)
                validQuantities = all(cellfun(@(q) any(strcmp(q,providedQuantities)),this.requiredQuantities));
            end


            if ~validQuantities
                msg2 = sprintf('Required quantities not provided by Dose Engine or selected Machine Dataset');
            end
            
            msg = strjoin({msg1,msg2},',');

            if isequal(msg(end),',')
                msg = msg(1:end-1);
            end

            avail = validRadMode && validQuantities;
        end
    end

    methods (Static)
        function classList = getAvailableModels(radiationMode,providedQuantities)
            %Get available models (optionally for radiationMode / machine combination)
            matRad_cfg = MatRad_Config.instance();

            %Use the root folder and the biomodel folder only
            folders = {fileparts(mfilename('fullpath'))};
            folders = [folders matRad_cfg.userfolders];

            persistent metaBioModels lastOptionalPaths
            if isempty(metaBioModels) || (~isempty(lastOptionalPaths) && ~isequal(lastOptionalPaths, folders))
                lastOptionalPaths = folders;
                metaBioModels = matRad_findSubclasses(meta.class.fromName(mfilename('class')),'folders',folders,'includeSubfolders',true);
            end
            
            classList = matRad_identifyClassesByConstantProperties(metaBioModels,'model','defaults',{'none'});

            if nargin > 0 
                callargs = {radiationMode};
                if nargin == 2 && iscellstr(providedQuantities)
                    callargs{end+1} = providedQuantities;
                end
                
                isAvailable = false(1,numel(classList));
                
                for i = 1:numel(classList)
                    modelInstance = classList(i).handle();
                    isAvailable(i) = all(modelInstance.isAvailable(callargs{:}));
                end

                classList = classList(isAvailable);
            end

            if isempty(classList)
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('No models found in paths %s',strjoin(folders,'\n'));
            end
        end
        
        function model = create(modelMetadata)
            if isa(modelMetadata,'matRad_BiologicalModel')
                model = modelMetadata;
                return;
            end

            matRad_cfg = MatRad_Config.instance();

            if ischar(modelMetadata) || isstring(modelMetadata)
                modelMetadata = struct('model',modelMetadata);
            end

            modelClassList = matRad_BiologicalModel.getAvailableModels();
            modelNames = {modelClassList.model};
            
            if ~isfield(modelMetadata,'model') || ~any(strcmp(modelNames,modelMetadata.model))
                matRad_cfg.dispWarning('Biological Model not found, creating Empty Model!');
                model = matRad_EmptyBiologicalModel();
                return;
            end

            usedModel = find(strcmp(modelNames,modelMetadata.model));
            
            if ~isscalar(usedModel)
                usedModel = usedModel(1);
            end
            
            modelClassInfo = modelClassList(usedModel);

            model = modelClassInfo.handle();

            modelMetadata = rmfield(modelMetadata,'model');
            
            %Now overwrite properties
            fields = fieldnames(modelMetadata);
            
            % iterate over all fieldnames and try to set the
            % corresponding properties inside the engine
            if matRad_cfg.isOctave
                c2sWarningState = warning('off','Octave:classdef-to-struct');                
            end
            
            for i = 1:length(fields)
                try
                    field = fields{i};
                    if matRad_ispropCompat(model,field)
                        model.(field) = matRad_recursiveFieldAssignment(model.(field),modelMetadata.(field),true);
                    else
                        matRad_cfg.dispWarning('Not able to assign property ''%s'' from bioModel struct to Biological Model!',field);
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
                            matRad_cfg.dispWarning('Not able to assign property from bioMOdel struct to Biological Model: %s',ME.message);
                        otherwise
                            matRad_cfg.dispWarning('Problem while setting up Biological Model from struct:%s %s',field,ME.message);
                    end
                end
            end
            
            if matRad_cfg.isOctave
                warning(c2sWarningState.state,'Octave:classdef-to-struct');                
            end
        end
    
        function model = validate(model,radiationMode,providedQuantities)
            %Make sure model is a validly created instance
            model = matRad_BiologicalModel.create(model);
            
            if nargin == 3
                [valid,msg] = model.isAvailable(radiationMode,providedQuantities);
            else
                [valid,msg] = model.isAvailable(radiationMode);
            end
            
            if ~valid
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Biological Model not valid: %s',msg);
            end
        end

        
    end

end