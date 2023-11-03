classdef matRad_RayTracer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties 
       lateralCutOff; 
    end
    
    
    methods
        function obj = matRad_RayTracer()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            matRad_cfg = MatRad_Config.instance();
            obj.lateralCutOff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;
        end
        
        function [alphas,l,rho,d12,ix] = traceRays(this,...
                isocenter, ...
                resolution, ...
                sourcePoints, ...
                targetPoints, ...
                cubes)
            
            %Default trivial implementation based on traceRay
            nRays = size(targetPoint,1);
            nSources = size(sourcePoint,1);
            
            if nSources ~= nRays && nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Number of source points (%d) needs to be one or equal to number of target points (%d)!',nSources,nRays);
            elseif nSources == 1
                sourcePoint = repmat(sourcePoint,nRays,1);
                nSources = nRays;
            end            
            
            for r = 1:nRays
                [alphas{r},l{r},rho{d12},ix{r}] = this.traceRay(isocenter,resolution,sourcePoints(r,:),targetPoints(r,:),cubes);
            end
            
            %pad with NaN values
            numval = cellfun(@numel,ix);
            maxnumval = max(numval);
            
            nanpad = @(x) [x(1:end), NaN(maxnumval - length(x),1)];
            
            alphas = cellfun(nanpad,alphas,'UniformOutput',false);
            l = cellfun(nanpad,l,'UniformOutput',false);
            ix = cellfun(nanpad,ix,'UniformOutput',false);
            
            for c = 1:numel(cubes)
                rho{c} = cellfun(nanpad,rho{c},'UniformOutput',false);
            end
            % Now the output should be consistent

        end
        
        function [alphas,l,rho,d12,ix] = traceRay(this,...
                isocenter, ...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes)
            error('Needs to be implemented!');   
        end
        
        function [radDepthsV,radDepthCube] = traceCube(this,stfElement,ct,V,rot_coordsV)
            matRad_cfg = MatRad_Config.instance();
            
            if ~isstruct(stfElement) || numel(stfElement) ~= 1
                matRad_cfg.dispError('The RayTracer does not accept stf struct arrays stf and only operates on a single field!');
            end
            
            %If no subset of voxels is specified, take all of them
            if nargin < 4
                V = transpose(1:prod(ct.cubeDim));
            end
                       
            %if we don't provide rotated patient coordinates we can compute
            %them here on our own
            if nargin < 5
                [yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,V);
                xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stfElement.isoCenter(1);
                yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stfElement.isoCenter(2);
                zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stfElement.isoCenter(3);
                coordsV        = [xCoordsV yCoordsV zCoordsV];
                
                % Get Rotation Matrix
                % Do not transpose matrix since we usage of row vectors &
                % transformation of the coordinate system need double transpose
                
                rotMat_system_T = matRad_getRotationMatrix(stfElement.gantryAngle,stfElement.couchAngle);
                
                % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
                rot_coordsV         = coordsV*rotMat_system_T;
                
                %
                rot_coordsV(:,1) = rot_coordsV(:,1)-stfElement.sourcePoint_bev(1);
                rot_coordsV(:,2) = rot_coordsV(:,2)-stfElement.sourcePoint_bev(2);
                rot_coordsV(:,3) = rot_coordsV(:,3)-stfElement.sourcePoint_bev(3);
            end
            
            % set up ray matrix direct behind last voxel
            rayMx_bev_y = max(rot_coordsV(:,2)) + max([ct.resolution.x ct.resolution.y ct.resolution.z]);
            rayMx_bev_y = rayMx_bev_y + stfElement.sourcePoint_bev(2);
            
            
            % set up list with bev coordinates for calculation of radiological depth
            coords = zeros(prod(ct.cubeDim),3);
            coords(V,:) = rot_coordsV;
            
            % calculate spacing of rays on ray matrix
            rayMxSpacing = 1/sqrt(2) * min([ct.resolution.x ct.resolution.y ct.resolution.z]);
            
            % define candidate ray matrix covering 1000x1000mm^2
            numOfCandidateRays = 2 * ceil(500/rayMxSpacing) + 1;
            candidateRayMx     = zeros(numOfCandidateRays);
            
            % define coordinates
            [candidateRaysCoords_X,candidateRaysCoords_Z] = meshgrid(rayMxSpacing*[floor(-500/rayMxSpacing):ceil(500/rayMxSpacing)]);
            
            % check which rays should be used
            for i = 1:stfElement.numOfRays
                
                ixCandidates = (candidateRaysCoords_X(:) - (1+rayMx_bev_y/stfElement.SAD)*stfElement.ray(i).rayPos_bev(1)).^2 + ...
                    (candidateRaysCoords_Z(:) - (1+rayMx_bev_y/stfElement.SAD)*stfElement.ray(i).rayPos_bev(3)).^2 ...
                    <= this.lateralCutOff^2;
                
                candidateRayMx(ixCandidates) = 1;
                
            end
            
            % set up ray matrix
            rayMx_bev = [candidateRaysCoords_X(logical(candidateRayMx(:))) ...
                rayMx_bev_y*ones(sum(candidateRayMx(:)),1) ...
                candidateRaysCoords_Z(logical(candidateRayMx(:)))];
            
            % Rotation matrix. Transposed because of row vectors
            rotMat_vectors_T = transpose(matRad_getRotationMatrix(stfElement.gantryAngle,stfElement.couchAngle));
            
            % rotate ray matrix from bev to world coordinates
            rayMx_world = rayMx_bev * rotMat_vectors_T;
            
            % criterium for ray selection
            raySelection = rayMxSpacing/2;
            
            %Trace all selected rays
            [~,l,rho,~,ixHitVoxel] = this.traceRays(stfElement.isoCenter, ...
                ct.resolution, ...
                stfElement.sourcePoint, ...
                rayMx_world, ...
                ct.cube);
            
            % find voxels for which we should remember this tracing because this is
            % the closest ray by projecting the voxel coordinates to the
            % intersection points with the ray matrix and checking if the distance
            % in x and z direction is smaller than the resolution of the ray matrix
            scale_factor = NaN(size(ixHitVoxel));
            valid_ix = ~isnan(ixHitVoxel);
            scale_factor(valid_ix) = (rayMx_bev_y - stfElement.sourcePoint_bev(2)) ./ ...
                   coords(ixHitVoxel(valid_ix),2);
            
            x_dist = NaN(size(ixHitVoxel));
            z_dist = NaN(size(ixHitVoxel));
            
            x_dist(valid_ix) = coords(ixHitVoxel(valid_ix),1).*scale_factor(valid_ix);
            x_dist = x_dist - rayMx_bev(:,1);
            
            z_dist(valid_ix) = coords(ixHitVoxel(valid_ix),3).*scale_factor(valid_ix);
            z_dist = z_dist - rayMx_bev(:,3);
            
            % Find indices
            ixRememberFromCurrTracing = x_dist > -raySelection & x_dist <= raySelection ...
                & z_dist > -raySelection & z_dist <= raySelection;
                       
            % set up rad depth cube for results
            radDepthCube = repmat({NaN(ct.cubeDim)},ct.numOfCtScen);
            radDepthV = cell(size(radDepthCube));
            
            for j = 1:ct.numOfCtScen
                rayDistances = l .* rho{j};
                rayWepl = cumsum(rayDistances,2) - rayDistances/2;
                
                radDepthCube{j}(ixHitVoxel(ixRememberFromCurrTracing)) = rayWepl(ixRememberFromCurrTracing);
                
                radDepthsV{j} = radDepthCube{j}(V);
            end
        end
        
        
    end
end

