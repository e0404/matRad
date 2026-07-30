classdef (Abstract) matRad_FluenceOptimizationFunction < matRad_OptimizationFunction
    % matRad_FluenceOptimizationFunction. Superclass for fluence objectives and
    % constraints.
    % This is the superclass for all objectives and constraints acting on the
    % bixel weight vector w itself instead of on a dose related quantity. No
    % dose is involved, so these functions take no part in the backprojection /
    % chain rule - their gradient is a gradient w.r.t. the weights already.
    %
    % A fluence function is tied to the beamlet layout, which it takes from the
    % dij in setupForDij. That is called once during the initialization of the
    % optimization problem (see matRad_getFluenceObjectives), not by the user.
    % isCompatible reports whether that has happened and whether the dimensions
    % still match.
    %
    % Note on what can be expressed: the dij only carries the *mapping* to the
    % stf (which bixel belongs to which beam), not the beamlet geometry itself.
    % Fluence functions built from a dij alone are therefore beam-wise. Anything
    % resolving the position of a bixel within its beam - e.g. smoothing that
    % couples neighbouring bixels - needs the stf and has to wait until the
    % optimization entry point receives it.
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

    properties (SetAccess = protected)
        numOfBixels = NaN       % number of bixels this function was set up for (NaN: not set up yet)
        bixelIx     = []        % global bixel indices the function acts on (empty: not set up yet)
    end

    methods

        function obj = matRad_FluenceOptimizationFunction(varargin)
            obj@matRad_OptimizationFunction(varargin{:});
        end

        function obj = setupForDij(obj, dij) %#ok<INUSD>
            % Build the dij dependent part of the function.
            % The default does nothing - functions that do not depend on the
            % beamlet layout are usable right after construction.
        end

        function [isCompat, msg] = isCompatible(obj, totalNumOfBixels)
            % Check whether the function can be evaluated for a weight vector
            % of the given length.

            msg = '';
            if isnan(obj.numOfBixels)
                isCompat = false;
                msg = sprintf('Fluence function ''%s'' has not been set up for a dij yet.', obj.name);
            elseif obj.numOfBixels ~= totalNumOfBixels
                isCompat = false;
                msg = sprintf(['Fluence function ''%s'' was set up for %d bixels, but the optimization ' ...
                               'problem has %d.'], obj.name, obj.numOfBixels, totalNumOfBixels);
            else
                isCompat = true;
            end
        end

    end

    methods (Static)

        function beams = getDefaultBeams(dij)
            % Beams a fluence function should act on by default. For a VMAT
            % plan only the FMO beams carry a fluence that is optimized; the
            % remaining control points are bounded to zero.
            if isfield(dij, 'isFMOBeam')
                beams = find(dij.isFMOBeam);
            else
                beams = 1:dij.numOfBeams;
            end
        end

        function beamBixelIx = getBeamBixels(dij, beams)
            % Global bixel indices of each requested beam, taken from the
            % bixel-to-beam mapping the dij carries. Works for any modality -
            % no assumption on bixels per ray.
            %
            % output:
            %   beamBixelIx:  cell array, one index vector per requested beam

            matRad_cfg = MatRad_Config.instance();

            if ~isfield(dij, 'beamNum') || numel(dij.beamNum) ~= dij.totalNumOfBixels
                matRad_cfg.dispError(['This dij does not carry a usable bixel-to-beam mapping (dij.beamNum). ' ...
                                      'Fluence optimization functions need it to identify the beams.']);
            end

            beamBixelIx = arrayfun(@(b) find(dij.beamNum(:) == b), beams(:)', 'UniformOutput', false);

            if any(cellfun(@isempty, beamBixelIx))
                matRad_cfg.dispWarning('Some of the requested beams carry no bixels and are ignored.');
                beamBixelIx = beamBixelIx(~cellfun(@isempty, beamBixelIx));
            end
        end

    end

end
