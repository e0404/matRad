classdef matRad_StfGeneratorBase < handle

    properties (Access = protected)
        machine
        pln
        cubeDim
        coordsY_vox 
        coordsX_vox 
        coordsZ_vox 
        voiTarget
        
    end

    properties
        visMode = true;
        %all propStf properties here except the brachy stf after the gens
        %are properly working
    end


    methods 
        function this = matRad_StfGeneratorBase(pln)
            matRad_cfg = MatRad_Config.instance();
            addpath(fullfile(matRad_cfg.matRadRoot));

            %after the gens are properly working implent a setDefaults from
            %pln function like in doseengines and assign all the propStf
            %properties here to override default attributes from pln
            
            this.pln = pln;

            
            
            if ~isfield(pln, 'propStf')
                matRad_cfg.dispError('no applicator information in pln struct');
            end
        end
    end

    methods 

        function stf = generate(this, ct, cst, visMode)  
            % Instance of MatRad_Config class
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('matRad: Generating stf struct... ');
            % load default parameters if not set
            this.pln = matRad_cfg.getDefaultProperties(this.pln, {'propOpt','propStf'});

            if nargin < 4
                this.visMode = 0;
            end

            % get machine
            try
                this.machine = matRad_loadMachine(this.pln);
            catch
                matRad_cfg.dispError('Could not find the following machine file: %s', fileName);
            end

            % Check Config
            if ~isfield(this.pln, 'multScen')
                matRad_cfg.dispWarning('No scenario model specified! Using nominal Scenario model!');
                this.pln.multScen = matRad_NominalScenario(ct);
            end

            this.initializePatientGeometry(ct, cst, visMode);
            stf = this.generateSourceGeometry(ct, cst, visMode);
        end
    end

    methods (Access = protected)

        function initializePatientGeometry(this, ct, cst, visMode)
            matRad_cfg = MatRad_Config.instance();
            
            % Initialize patient geometry
            V = [];
            ct = matRad_calcWaterEqD(ct,this.pln);  

            isTarget = cellfun(@(voiType) isequal(voiType, 'TARGET'), cst(:,3));
            if ~any(isTarget)
                matRad_cfg.dispError('No target found in cst. Please designate at least one VOI as ''TARGET''!');
            end

            hasObjective = ~cellfun(@isempty, cst(:,6));
            useTargetForBixelPlacement = isTarget & hasObjective;

            if ~any(useTargetForBixelPlacement)
                matRad_cfg.dispWarning('No Objectives / Constraints assigned to targets. All targets will be considered for Bixel placement!');
                useTargetForBixelPlacement(isTarget) = true;
            end

            % Now add all used target voxels to the voxel list
            for i = 1:size(cst, 1)
                if useTargetForBixelPlacement(i)
                    V = [V; cst{i,4}{1}];  
                end
            end

            % Remove double voxels
            V = unique(V);
            % generate voi cube for targets
            this.voiTarget = zeros(ct.cubeDim);
            this.voiTarget(V) = 1;

            % add margin information
            addmarginBool = matRad_cfg.defaults.propStf.addMargin;
            if isfield(this.pln, 'propStf') && isfield(this.pln.propStf, 'addMargin')
                addmarginBool = this.pln.propStf.addMargin;
            end

            % Margin info
            if addmarginBool
                pbMargin = this.getPbMargin();

                % Assumption for range uncertainty
                assumeRangeMargin = this.pln.multScen.maxAbsRangeShift + this.pln.multScen.maxRelRangeShift + pbMargin;

                % add margin - account for voxel resolution, the maximum shift scenario and the current bixel width.
                margin.x = max([ct.resolution.x max(abs(this.pln.multScen.isoShift(:,1)) + assumeRangeMargin)]);
                margin.y = max([ct.resolution.y max(abs(this.pln.multScen.isoShift(:,2)) + assumeRangeMargin)]);
                margin.z = max([ct.resolution.z max(abs(this.pln.multScen.isoShift(:,3)) + assumeRangeMargin)]);

                this.voiTarget = matRad_addMargin(this.voiTarget, cst, ct.resolution, margin, true);
                V = find(this.voiTarget > 0);
            end

            % throw error message if no target is found
            if isempty(V)
                matRad_cfg.dispError('Could not find target.');
            end

            % Convert linear indices to 3D voxel coordinates
            [this.coordsY_vox, this.coordsX_vox, this.coordsZ_vox] = ind2sub(ct.cubeDim, V);

            if isempty(this.coordsX_vox) || isempty(this.coordsY_vox) || isempty(this.coordsZ_vox)
                matRad_cfg.dispWarning('coordsXYZ are empty, boundary cannot be computed.');   % they AREN'T EMPTY here problem is in brachyStfGen
            end

            % take only voxels inside patient
            V = [cst{:,4}];
            V = unique(vertcat(V{:}));

            % ignore densities outside of contours
            eraseCtDensMask = ones(prod(ct.cubeDim), 1);
            eraseCtDensMask(V) = 0;
            for i = 1:ct.numOfCtScen
                ct.cube{i}(eraseCtDensMask == 1) = 0;
            end

            
        end

        function pbMargin = getPbMargin(this)
            pbMargin = 0;
        end
    end

    methods (Access = protected)
        % the actual calculation method which returns the final stf struct.
        % Needs to be implemented in non-abstract subclasses.
        % (Internal logic is often split into multiple methods in order to
        % make the whole calculation more modular)
        function stf = generateSourceGeometry(this, ct, cst, visMode)
            throw(MException('MATLAB:class:AbstractMember','Abstract function generateSourceGeometry of your StfGenerator needs to be implemented!'));
        end
    end
end
