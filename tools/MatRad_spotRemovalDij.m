classdef MatRad_spotRemovalDij < handle
    % MatRad_spotRemovalDij class definition
    %
    %
    % References
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties

        matRad_cfg = MatRad_Config.instance();

        removalMode = 'relative'; % 'relative' is the only mode for now
        propSpotRemoval;

        dij;
        cst;
        pln;
        stf;

        weights;
        newSpots;
        newWeights;

        numOfRemovedSpots;

    end

    methods
        function obj = MatRad_spotRemovalDij(dij,w)

            obj.reset();

            % I don't know if this is good practice to copy a very large struct
            obj.dij = dij;
            obj.weights = w;

        end


        function reset(obj)

            %Set all default properties for spot removal
            obj.setDefaultProperties();

        end


        function obj = setDefaultProperties(obj)

            obj.propSpotRemoval.relativeThreshold = 0.03;
            obj.propSpotRemoval.absoluteThreshold = 2;

        end


        function resultGUI = reoptimize(obj,cst,pln)

            if isempty(obj.cst) || exist('cst','var')
                obj.cst = cst;
            end
            if isempty(obj.pln) || exist('pln','var')
                obj.pln = pln;
            end

            if isempty(obj.newWeights)
                obj.calcNewSpots();
            end

            resultGUI = matRad_fluenceOptimization(obj.getDij,obj.cst,obj.pln,obj.newWeights);

        end


        function obj = calcNewSpots(obj)

            switch obj.removalMode
                case 'relative'
                    obj.newSpots = obj.weights > obj.propSpotRemoval.relativeThreshold * mean(obj.weights);
                    obj.matRad_cfg.dispInfo([num2str(sum(~obj.newSpots)),'/',num2str(numel(obj.newSpots)) ,' spots have been removed below ',num2str(100*obj.propSpotRemoval.relativeThreshold),'% of the mean weight.\n'])
                case 'absolute'
                    obj.newSpots = obj.weights > obj.propSpotRemoval.absoluteThreshold;
                    obj.matRad_cfg.dispInfo([num2str(sum(~obj.newSpots)),'/',num2str(numel(obj.newSpots)) ,' spots have been removed below thres=',num2str(obj.propSpotRemoval.absoluteThreshold),'.\n'])
                otherwise
                    obj.matRad_cfg.dispWarning(['Removal mode ' obj.removalMode ' not implemented, no spots have been removed.']);
            end

            obj.numOfRemovedSpots = sum(~obj.newSpots);
            obj.newWeights = obj.weights(obj.newSpots);

        end


        function stf = getStf(obj,stf)

            if isempty(obj.stf) || exist('stf','var')
                obj.stf = stf;
            end

            if ~isempty(obj.newWeights) && ~isempty(obj.stf)
                stf = obj.stf;
                [~,beamNumIdx] = unique(obj.dij.beamNum);
                beamNumIdx = [0;beamNumIdx(2:end)-1;obj.dij.totalNumOfBixels];

                for b = 1:obj.dij.numOfBeams
                    currRaysInBeam = obj.dij.rayNum(beamNumIdx(b)+1:beamNumIdx(b+1));
                    currBixelsInRay = obj.dij.bixelNum(beamNumIdx(b)+1:beamNumIdx(b+1));
                    [rayCnt,rayIdx] = unique(currRaysInBeam);

                    numOfBixelsPerRay = groupcounts(currRaysInBeam);
                    cutRays = ismember([1:obj.dij.numOfRaysPerBeam(b)]',rayCnt);
                    if any(~cutRays)
                        stf(b).ray = stf(b).ray(cutRays);
                        stf(b).numOfRays = sum(cutRays);
                    end
                    bixelCurrRay = cell(1,stf(b).numOfRays);
                    for i = 1:stf(b).numOfRays
                        bixelCurrRay{i} = currBixelsInRay(rayIdx(i):rayIdx(i)+numOfBixelsPerRay(i)-1);
                    end
                    for f = 1:stf(b).numOfRays
                        stf(b).ray(f).energy = stf(b).ray(f).energy(bixelCurrRay{f});
                        stf(b).ray(f).focusIx = stf(b).ray(f).focusIx(bixelCurrRay{f});
                        stf(b).ray(f).rangeShifter = stf(b).ray(f).rangeShifter(bixelCurrRay{f});
                    end
                    stf(b).numOfBixelsPerRay = numOfBixelsPerRay';
                    stf(b).totalNumOfBixels = sum(stf(b).numOfBixelsPerRay);
                end
            end

        end


        function dij = getDij(obj)

            if isempty(obj.newWeights)
                obj.calcNewSpots();
            end

            dij = obj.dij;
            if ~isempty(obj.newWeights)
                dij.cutWeights = obj.newWeights;

                dij.bixelNum = dij.bixelNum(obj.newSpots);
                dij.rayNum = dij.rayNum(obj.newSpots);
                dij.beamNum = dij.beamNum(obj.newSpots);
                dij.totalNumOfBixels = sum(obj.newSpots);

                dij.physicalDose{1} = dij.physicalDose{1}(:,obj.newSpots);
                if isfield(dij,'mAlphaDose')
                    dij.mAlphaDose{1} = dij.mAlphaDose{1}(:,obj.newSpots);
                    dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}(:,obj.newSpots);
                end
                if isfield(dij,'mLETDose')
                    dij.mLETDose{1} = dij.mLETDose{1}(:,obj.newSpots);
                end
                [~,beamNumIdx] = unique(dij.beamNum);
                beamNumIdx = [0;beamNumIdx(2:end)-1;dij.totalNumOfBixels];

                for b = 1:dij.numOfBeams
                    currRaysInBeam = dij.rayNum(beamNumIdx(b)+1:beamNumIdx(b+1));
                    [rayCnt,~] = unique(currRaysInBeam);

                    dij.numOfRaysPerBeam(b) = numel(rayCnt);
                end

                dij.totalNumOfRays = sum(dij.numOfRaysPerBeam);
                dij.numOfRemovedSpots = sum(~obj.newSpots);
            end

        end


        function weights = getWeights(obj)

            if isempty(obj.newWeights)
                obj.calcNewSpots();
            end
            if ~isempty(obj.newWeights)
                weights = obj.newWeights;
            end

        end


        function weightsLogic = getLogical(obj)

            if isempty(obj.newWeights)
                obj.calcNewSpots();
            end
            if ~isempty(obj.newWeights)
                weightsLogic = obj.newSpots;
            end

        end

    end
end


