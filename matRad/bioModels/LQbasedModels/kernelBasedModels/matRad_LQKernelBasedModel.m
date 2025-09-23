classdef (Abstract) matRad_LQKernelBasedModel < matRad_LQBasedModel
% This is an Abstract class implementing any Kernel-Based model.
% This class of models rely on the interpolation of the variable alpha/beta
% kernels presen in the base data.
% The different tissue parameters are specified by the cst structure as
% tissue classes.
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

    properties (Abstract, Constant)
        kernelQuantities;
    end

    methods
        function this = matRad_LQKernelBasedModel()
            this@matRad_LQBasedModel();
        end
    end

    methods %(Static)
        
        function vTissueIndex = getTissueInformation(this,machine, cstDownsampled, dij,vAlphaX, ~, VdoseGrid, VdoseGridScenIx)

            matRad_cfg = MatRad_Config.instance();

            numOfCtScen = numel(vAlphaX);
            
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

    methods (Static)


        function [alphaX, betaX] = getAvailableTissueParameters(pln)
            
            % load machine
            machine = matRad_loadMachine(pln);
            if isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
                alphaX = machine.data(1).alphaX;
                betaX  = machine.data(1).betaX; 
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('The selected biological model requires AlphaX and BetaX to be set in the machine file but none was found.');
            end

        end
    end
end