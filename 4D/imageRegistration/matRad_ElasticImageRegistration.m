classdef  matRad_ElasticImageRegistration < matRad_ImageRegistration
    
    % matRad elastic calculation
    %
    % call
    %   obj = matRad_elasticImageRegistration(ct,cst,refScen)
    %   obj = matRad_elasticImageRegistration(ct,cst,refScen,metadata)
    %
    % input
    %   ct:             matRad ct struct
    %   cst:            matRad cst struct
    %   metadata:       struct of metadata
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2022 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    properties (Constant)
        name = 'Elastic Registration';
    end
    
    properties
        ct
        cst
        refScen
        metadata
    end
    
    methods (Access = public)
        
        
        function obj = matRad_ElasticImageRegistration(ct,cst,refScen,metadata)
            
            if ~isfield(metadata,'dvfType')
                metadata.dvfType = 'pull';
            end
            
            if ~isfield(metadata,'nItera')
                metadata.nItera = 100; %Default number of iterations
            end
            
            if ~isfield(metadata,'pyramLevels')
                metadata.pyramLevels = 1; %Default number of multi-resolution image pyramid levels
            end
            
            if ~isfield(metadata,'smoothLevels')
                metadata.smoothLevels = 1; %Default smooth levels
            end
            
            if ~exist('refScen','var') || isempty(refScen)
                refScen = 1;
            end
            
            obj.ct = ct;
            obj.cst = cst;
            obj.refScen = refScen;
            obj.metadata=metadata;
            
        end
        
        % Calculate the deformation vector fields
        function [ct,cst] = calcDVF(obj)
            
            % Non rigid registration demons-based. Calculates the DVF(Displacement
            % Vector Field) that models the transformation.
            
            nScen = obj.ct.numOfCtScen;
            obj.ct.dvf = cell(1,nScen);
            obj.ct.dvfType = obj.metadata.dvfType;
            obj.ct.refScen = obj.refScen;
            
            switch obj.metadata.dvfType
                case 'pull'
                    for scen = 1:nScen
                        fprintf('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(obj.ct.cubeHU{obj.refScen},obj.ct.cubeHU{scen},obj.metadata.nItera,'PyramidLevels',obj.metadata.pyramLevels,'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
                case 'push'
                    for scen = 1:nScen
                        fprintf('Registering scenario %d.\n',scen);
                        [obj.ct.dvf{scen},~] = imregdemons(obj.ct.cubeHU{scen},obj.ct.cubeHU{obj.refScen},obj.metadata.nItera,'PyramidLevels',obj.metadata.pyramLevels,'AccumulatedFieldSmoothing',obj.metadata.smoothLevels);
                        obj.ct.dvf{scen} = permute(obj.ct.dvf{scen},[4 1 2 3]);
                    end
            end
            
            ct=obj.ct;
            cst=obj.cst;
            
        end
        
        % Propagate contourns
        function [ct,cst] = propContours(obj)
            
            if ~strcmp(obj.metadata.dvfType,'push')
                error('error propagation requires push dvfs');
            end
            
            if ~isfield(obj.ct,'dvf') || ~isfield(obj.ct,'dvfType')
                [obj.ct,obj.cst]=calcDVF(obj);
            end
            
            [numOfStruct, ~] = size(obj.cst);
            for structure = 1:numOfStruct
                
                if ~isempty(obj.cst{structure,4}{1,1})
                    
                    % Obtaining the fixed cubic structure from the linear indices
                    cubeHU_fixed = zeros(obj.ct.cubeDim);
                    cst_fixed = obj.cst{structure,4}{1,1};
                    [x,y,z] = ind2sub(obj.ct.cubeDim,cst_fixed);
                    
                    % The HU value of the corresponding tomography is assigned to each position
                    for j=1:length(x)
                        cubeHU_fixed(x(j),y(j),z(j)) = 1; %obj.cubeHU{1}(x(j),y(j),z(j));
                    end
                    
                    % The DVF transformation is applied and the linear values are found
                    fprintf('Propagating countorns of structure %d \n',structure);
                    for scen = 2:obj.ct.numOfCtScen
                        cubeHU_estimated = imwarp(cubeHU_fixed, permute(obj.ct.dvf{scen_moving},[2 3 4 1]));
                        obj.cst{structure,4}{1,scen} = find(cubeHU_estimated);
                    end
                else
                    fprintf('Deleting empty structure %d \n',structure);
                    obj.cst(structure,:) = [];
                end
            end
            
            ct=obj.ct;
            cst=obj.cst;
            
        end
        
    end
    
end

