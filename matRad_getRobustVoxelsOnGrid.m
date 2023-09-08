function [includeMask] = matRad_getRobustVoxelsOnGrid(cst, doseGrid, selectionMode)
% matRad function to get mask of the voxels (on dose grid) that are
% included in structures used for robust optimization. Dose calculation in
% the non nonminal scenarios will anly be performed on the voxels sset to 1
% in the mask 
%
% call
%   includeMask = matRad_getRobustVoxelsOnGrid(cst,doseGrid,VdoseGrid,selectionMode)
%
% input
%   cst:                    cst (voxel indexes included in cst{:,4} are referred to a cube of dimensions doseGrid.dimensions)
%   doseGrid:               doseGrid struct containing field doseGrid.dimensions
%   selectionMode:          define wich method to apply to select the cst
%                           structures to include. Choices are: 
%                               none            does nothing, all voxels will be included                                                                                  
%                               all             includes all the structures in the cst%                               
%                               targetOnly      only includes the voxels in structures labeld as target                      
%                               oarsOnly        only includes the voxels in structures labeld as oars    
%                               objectivesOnly  only includes the voxels in structures with at least one objective/constraint    
%                               TODO:           all structurtes with robustness objectives
%                               
% 
% output
% 
%   includeMask:            logical array #voxels in dose grid
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    matRad_cfg = MatRad_Config.instance();
    
    matRad_cfg.dispInfo('Disabling dose calculation in voxels outside of ROIs in robustness scenarios');
    cst  = matRad_setOverlapPriorities(cst);

    selectedCstStructs = [];
    
    includeMask = cell(size(cst{1,4},2),1);
    includeMask(:) = {zeros(prod(doseGrid.dimensions),1)};
        
    if isequal(selectionMode , 'none')
        for ctScenIdx=1:size(includeMask,2)
            includeMask{ctScenIdx}(:) = 1;
        end
    
    else
        %structures that are selected here will be included in dose calculation
        %over the robust scenarios

        switch selectionMode
        
            case 'all'
                selectedCstStructs = [1:size(cst,1)];
            case 'targetOnly'
                selectedCstStructs = find(cellfun(@(x) strcmp(x,'TARGET'), [cst(:,3)]));
            case 'objectivesOnly'
                for i=1:size(cst,1)
                    if numel(cst{i,6})>0
                        selectedCstStructs = [selectedCstStructs, i];
                    end
                end
            case 'oarsOnly'
                selectedCstStructs = find(cellfun(@(x) strcmp(x,'OAR'), [cst(:,3)]));    
            otherwise
        
        end
    
        
        %loop over all cst sturctures 
        for i=1:size(cst,1)
        
            if ~isempty(cst{i,4}{1})
                
                if numel(cst{i,6}) > 0
                    %loop over obj/constraint functions
                    for j = 1:numel(cst{i,6})
        
    
                        obj = cst{i,6}{j};
        
                        if ~isa(obj,'matRad_DoseOptimizationFunction')
                            try
                                obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                            catch
                                matRad_cfg.dispError('cst{%d,6}{%d} is not a valid Objective/constraint! Remove or Replace and try again!',i,j);
                            end
                        end
    
                        robustness = obj.robustness;
    
                        if any(intersect(i, selectedCstStructs))
                            for ctIdx=1:size(cst{i,4},2)
                                includeMask{ctIdx}(cst{i,4}{ctIdx}) = 1;
                            end
    
                            if isequal(robustness, 'none')
                                matRad_cfg.dispWarning('Dose calculation is performed on cst structure %s, but this structure has no robustness.', cst{i,2});
                            end
                        else
                            if ~isequal(robustness, 'none')
                                matRad_cfg.dispError('Trying to clear voxels in cst structure: %s, but this structure has robustness requirements', cst{i,2});
                            end
                        end
                    end
    
                else %numel(cst{i,6}) < 0
                    if any(intersect(i, selectedCstStructs))
                        matRad_cfg.dispWarning('Trying to clear voxels in cst structure: %s, but this structure does not have any objective or constraint', cst{i,2}');
                    end
                end %numel(cst{i,6}) > 0
            end %if cst{i,4} not empty
        

        end %for loop over cst
    end

    for ctIdx=1:size(cst{1,4},2)
        includeMask{ctIdx} = includeMask{ctIdx};
    end
    
    matRad_cfg.dispInfo('\n done\n');
end