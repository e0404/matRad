classdef matRad_FluenceVariance < FluenceObjectives.matRad_FluenceObjective
    % matRad_FluenceVariance: penalizes the fluence variance within each beam
    %
    %   Beam-wise modulation reduction. With m_b the mean weight of beam b and
    %   N the number of bixels the objective acts on,
    %
    %       f(w) = 1/N * sum_b sum_{i in b} (w_i - m_b)^2  /  normalization
    %
    %   which is the quadratic form w'Sw with S block diagonal and
    %   S_bb = (I - 11'/n_b)/N - see the static getQuadraticForm. A flat
    %   fluence per beam is free, any modulation costs.
    %
    %   Bixels of different beams are never coupled, and no assumption is made
    %   about where a bixel sits inside its beam: the objective only needs the
    %   bixel-to-beam mapping the dij carries, so it works for photons and
    %   particles alike. Smoothing that couples *neighbouring* bixels needs the
    %   beamlet geometry from the stf and is not available here - see
    %   matRad_FluenceOptimizationFunction.
    %
    %   Parameters:
    %     normalization  'relative' (default) each beam's variance is divided
    %                    by that beam's OWN mean square and the beams are
    %                    averaged, i.e. the mean squared coefficient of
    %                    variation. Dimensionless, invariant under a rescaling
    %                    of the fluence, so penalties transfer between plans.
    %                    'absolute'  the plain quadratic form; its value scales
    %                    with the square of the fluence magnitude
    %
    %   The per-beam denominator matters: normalizing by a mean square pooled
    %   over all beams would put the differences between the beam levels into
    %   the denominator, since
    %       sum_i w_i^2 = sum_b sum_{i in b} (w_i - m_b)^2 + sum_b n_b m_b^2 ,
    %   so the objective could be lowered by driving the beams apart in
    %   intensity instead of by smoothing. On an arc that concentrates the
    %   plan on a few control points, which is the opposite of what this
    %   objective is for.
    %
    % References
    %   Quadratic fluence regularization of this type is standard in fluence
    %   map optimization, e.g. https://doi.org/10.1088/0031-9155/50/8/010
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2026 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Fluence Variance'
        parameterNames = {'normalization'}
        parameterTypes = {{'relative', 'absolute'}}
    end

    properties (Access = public)
        parameters = {1}
        penalty = 1
        beams = []          % dij beam indices to act on; empty selects the default (FMO) beams
    end

    properties (SetAccess = protected)
        beamBixelIx = {}    % bixel indices per selected beam, filled by setupForDij
    end

    methods

        function obj = matRad_FluenceVariance(penalty, normalization)
            % Construct the objective. The beam partition is filled in later
            % from the dij, see setupForDij.

            if nargin == 1 && isstruct(penalty)
                inputStruct = penalty;
                initFromStruct = true;
            else
                initFromStruct = false;
                inputStruct = [];
            end

            obj@FluenceObjectives.matRad_FluenceObjective(inputStruct);

            if ~initFromStruct
                if nargin >= 1 && isscalar(penalty)
                    obj.penalty = penalty;
                end
                if nargin >= 2 && ~isempty(normalization)
                    obj.parameters{1} = obj.matchOption(normalization, obj.parameterTypes{1}, 'normalization');
                end
            end
        end

        function normalization = getNormalization(obj)
            normalization = obj.parameterTypes{1}{obj.parameters{1}};
        end

        function obj = setupForDij(obj, dij)
            % Partition the bixels by beam. Called during the initialization of
            % the optimization problem.

            matRad_cfg = MatRad_Config.instance();

            usedBeams = obj.beams;
            if isempty(usedBeams)
                usedBeams = matRad_FluenceOptimizationFunction.getDefaultBeams(dij);
            end

            obj.beamBixelIx = matRad_FluenceOptimizationFunction.getBeamBixels(dij, usedBeams);
            obj.bixelIx     = sort(vertcat(obj.beamBixelIx{:}));
            obj.numOfBixels = dij.totalNumOfBixels;

            matRad_cfg.dispInfo('Fluence variance objective: %s, %d of %d bixels over %d beam(s)\n', ...
                                obj.getNormalization(), numel(obj.bixelIx), obj.numOfBixels, numel(obj.beamBixelIx));
        end

        function fFluence = computeFluenceObjectiveFunction(obj, w)
            fFluence = obj.evaluate(w);
        end

        function fFluenceGrad = computeFluenceObjectiveGradient(obj, w)
            % the gradient is always returned as a column, matching the weight
            % gradient of the optimization problem
            [~, fFluenceGrad] = obj.evaluate(w);
        end

        function s = struct(obj)
            s = struct@FluenceObjectives.matRad_FluenceObjective(obj);
            s.beams = obj.beams;
        end

    end

    methods (Access = protected)

        function [fVal, fGrad] = evaluate(obj, w)
            % Shared evaluation of the objective and its gradient, beam by
            % beam. O(number of bixels) - the quadratic form is never
            % materialized.

            matRad_cfg = MatRad_Config.instance();

            if isempty(obj.beamBixelIx)
                matRad_cfg.dispError('Fluence objective ''%s'' has not been set up for a dij!', obj.name);
            end

            switch obj.getNormalization()
                case {'relative', 'absolute'}
                    isRelative = strcmp(obj.getNormalization(), 'relative');
                otherwise
                    matRad_cfg.dispError('Invalid normalization ''%s''!', obj.getNormalization());
            end

            w = w(:);
            numActive = numel(obj.bixelIx);
            numBeams  = numel(obj.beamBixelIx);

            fVal  = 0;
            fGrad = zeros(numel(w), 1);

            for i = 1:numBeams
                ix = obj.beamBixelIx{i};
                wBeam = w(ix);
                numInBeam = numel(wBeam);

                deviation = wBeam - mean(wBeam);
                % the mean correction cancels in the gradient of the variance,
                % since sum_j d_j = 0
                beamVar     = (deviation' * deviation) / numInBeam;
                beamVarGrad = 2 * deviation / numInBeam;

                if isRelative
                    % normalize with the beam's OWN mean square. A denominator
                    % pooled over all beams would also contain the differences
                    % between the beam levels, and could then be inflated by
                    % driving the beams apart in intensity - which lowers the
                    % objective without smoothing anything.
                    beamMeanSq = (wBeam' * wBeam) / numInBeam;
                    if beamMeanSq <= 0
                        % a beam without fluence carries no modulation to speak
                        % of; keep the objective (and its gradient) finite
                        continue
                    end

                    fVal = fVal + beamVar / beamMeanSq / numBeams;
                    fGrad(ix) = fGrad(ix) + (beamVarGrad / beamMeanSq - ...
                                             beamVar / beamMeanSq^2 * (2 * wBeam / numInBeam)) / numBeams;
                else
                    % plain pooled quadratic form, see getQuadraticForm
                    fVal = fVal + beamVar * numInBeam / numActive;
                    fGrad(ix) = fGrad(ix) + beamVarGrad * numInBeam / numActive;
                end
            end
        end

    end

    methods (Access = private)

        function ix = matchOption(obj, value, options, what) %#ok<INUSL>
            matRad_cfg = MatRad_Config.instance();
            if isnumeric(value) && isscalar(value) && value >= 1 && value <= numel(options)
                ix = value;
                return
            end
            ix = find(strcmpi(value, options), 1);
            if isempty(ix)
                matRad_cfg.dispError('Invalid %s ''%s''! Valid values are: %s', ...
                                     what, char(value), strjoin(options, ', '));
            end
        end

    end

    methods (Static)

        function [quadForm, activeBixelIx] = getQuadraticForm(dij, beams)
            % Explicit quadratic form S with f = w'Sw for the 'absolute'
            % normalization, without needing an objective instance. The
            % 'relative' normalization divides each beam block by that beam's
            % own mean square and is therefore not a single quadratic form.
            %
            % The objective itself never builds this - it evaluates the same
            % expression in O(numOfBixels). S is provided for inspection and
            % for reuse elsewhere, and is dense within each beam block, so it
            % is only practical for moderate beam sizes.
            %
            % call:
            %   S = FluenceObjectives.matRad_FluenceVariance.getQuadraticForm(dij)
            %   S = FluenceObjectives.matRad_FluenceVariance.getQuadraticForm(dij,beams)

            matRad_cfg = MatRad_Config.instance();

            if nargin < 2 || isempty(beams)
                beams = matRad_FluenceOptimizationFunction.getDefaultBeams(dij);
            end

            beamBixelIx = matRad_FluenceOptimizationFunction.getBeamBixels(dij, beams);
            activeBixelIx = sort(vertcat(beamBixelIx{:}));
            numActive = numel(activeBixelIx);

            numEntries = sum(cellfun(@(ix) numel(ix)^2, beamBixelIx));
            if numEntries > 1e7
                matRad_cfg.dispError(['The explicit quadratic form would need %.3g entries - it is dense per ' ...
                                      'beam. Use the objective, which evaluates it without building S.'], numEntries);
            end

            rowIx = zeros(numEntries, 1);
            colIx = zeros(numEntries, 1);
            entries = zeros(numEntries, 1);
            offset = 0;

            for i = 1:numel(beamBixelIx)
                ix = beamBixelIx{i};
                n = numel(ix);
                block = (eye(n) - ones(n) / n) / numActive;

                [rr, cc] = ndgrid(ix, ix);
                rowIx(offset + (1:n^2)) = rr(:);
                colIx(offset + (1:n^2)) = cc(:);
                entries(offset + (1:n^2)) = block(:);
                offset = offset + n^2;
            end

            quadForm = sparse(rowIx, colIx, entries, dij.totalNumOfBixels, dij.totalNumOfBixels);
        end

    end

end
