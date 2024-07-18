classdef matRad_LQKernelBasedModel < matRad_LQBasedModel

    properties

    end

    methods
        function this = matRad_LQKernelBasedModel()
            this@matRad_LQBasedModel();

        end
    end

    methods %(Static)
        
        function vTissueIndex = getTissueInformation(this,machine, cst, dij,vAlphaX, ~, VdoseGrid, VdoseGridScenIx)

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

            if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
                for i = 1:size(cstDownsampled,1)
                    % check if cst is compatiable
                    if ~isempty(cstDownsampled{i,5}) && isfield(cstDownsampled{i,5},'alphaX') && isfield(cstDownsampled{i,5},'betaX')
     
                        % check if base data contains alphaX and betaX
                        % (should already be checked before to instantiate the model class)
                        IdxTissue = find(ismember(machine.data(1).alphaX,cstDownsampled{i,5}.alphaX) & ...
                            ismember(machine.data(1).betaX,cstDownsampled{i,5}.betaX));
     
                        % check consitency of biological baseData and cst settings
                        if ~isempty(IdxTissue)

                            for s = 1:numOfCtScen
                                tmpScenVdoseGrid = VdoseGrid(VdoseGridScenIx{s});
                                isInVdoseGrid = ismember(tmpScenVdoseGrid,cstDownsampled{i,4}{s});
                                vTissueIndex{s}(isInVdoseGrid) = IdxTissue;
                            end
                        else
                            matRad_cfg.dispError('Biological base data and cst are inconsistent!');
                        end
     
                    else
                        for s = 1:numOfCtScen
                            vTissueIndex{s}(:) = 1;
                        end
                        matRad_cfg.dispWarning('\tTissue type of %s was set to 1\n',cstDownsampled{i,2});
                    end
                end
                
                matRad_cfg.dispInfo('done.\n');
             else
                matRad_cfg.dispError('Base data is missing alphaX and/or betaX!');
            end

        end
    end
end