classdef matRad_RayTracerSiddon < matRad_RayTracer
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (SetAccess = private)
        sourcePoint;
        targetPoint;
        rayVec;
        
        xPlanes;
        yPlanes;
        zPlanes;
        numPlanes;
        resolution;
        cubeDim;
        
    end
    
    methods
        function this = matRad_RayTracerSiddon()
            
        end
        
        function [alphas,l,rho,d12,ix] = traceRay(this,...
                isocenter, ...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes)
            
            %traceRay Traces an individual ray
            %   Detailed explanation goes here            
            
            nRays = size(targetPoint,1);
            nSources = size(sourcePoint,1);
            
            if nRays ~= 1 || nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Number of target Points and source points needs to be equal to one! If you want to trace multiple rays at once, use traceRays instead!',nSources,nRays);
            end
            
            [alphas,l,rho,d12,ix] = this.traceRays(isocenter,resolution,sourcePoint,targetPoint,cubes);             
        end        
        
        function [alphas,l,rho,d12,ix] = traceRays(this,...
                isocenter, ...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes)
            
            nRays = size(targetPoint,1);
            nSources = size(sourcePoint,1);
            
            if nSources ~= nRays && nSources ~= 1
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Number of source points (%d) needs to be one or equal to number of target points (%d)!',nSources,nRays);
            elseif nSources == 1
                sourcePoint = repmat(sourcePoint,nRays,1);
                nSources = nRays;
            end
            
            this.sourcePoint = sourcePoint + isocenter;
            this.targetPoint = targetPoint + isocenter;
            
            this.initializeGeometry(resolution,size(cubes{1}));

            
            
            % eq 7 & 8
            % Calculate relative distances (alphas) at which intersections
            % occur
            alphas = this.computeAllAlphas();
                   
            % eq 11
            % Calculate the distance from source to target point.
            d12 = vecnorm(this.rayVec,2,2);
            
            
            % eq 10
            % Calculate the voxel intersection length.
            tmpDiff = diff(alphas,1,2);
            
            l = d12.*tmpDiff;
            alphas_mid = alphas(:,1:end-1) + 0.5*tmpDiff;
            
            
            
            % eq 12
            % Calculate the voxel indices: first convert to physical coords
            % and convert to voxel indices
            i = round((this.sourcePoint(:,1) + alphas_mid.*this.rayVec(:,1))./resolution.x);
            j = round((this.sourcePoint(:,2) + alphas_mid.*this.rayVec(:,2))./resolution.y);
            k = round((this.sourcePoint(:,3) + alphas_mid.*this.rayVec(:,3))./resolution.z);
            
            
            % Handle numerical instabilities at the borders.
            i(i<1) = 1;   i(i>length(this.xPlanes)-1) = length(this.xPlanes) - 1;
            j(j<1) = 1;   j(j>length(this.yPlanes)-1) = length(this.yPlanes) - 1;
            k(k<1) = 1;   k(k>length(this.zPlanes)-1) = length(this.zPlanes) - 1;
            
            valIx = ~isnan(alphas_mid);
            
            %In Matlab direct assignment with sub2ind would work, Octave 
            %however does not like NaN values in the subscripts
            ix = NaN(size(valIx));
            ix(valIx) = sub2ind([this.numPlanes(2), this.numPlanes(1), this.numPlanes(3)]-1,j(valIx),i(valIx),k(valIx));            
            
            for i = 1:numel(cubes)
                rho{i} = NaN(size(valIx));
                rho{i}(valIx) = cubes{i}(ix(valIx));
            end
            
            
        end
        

    end
    
    methods (Access = protected)       
        function alphas = computeAllAlphas(this)

            % Here we setup grids to enable logical indexing when computing
            % the alphas along each dimension. All alphas between the
            % minimum and maximum index will be computed, with additional
            % exclusion of singular plane occurences (max == min)
            % All values out of scope will be set to NaN.
            
            nRays = size(this.rayVec,1);
            
            % eq 4
            % Calculate parametrics values of \alpha_{min} and \alpha_{max} for every
            % axis, intersecting the ray with the sides of the CT.
            aX_1 = (this.xPlanes(1) - this.sourcePoint(:,1)) ./ this.rayVec(:,1);
            aX_end = (this.xPlanes(end) - this.sourcePoint(:,1)) ./ this.rayVec(:,1);
            
            tmpIx = this.targetPoint(:,1) == this.sourcePoint(:,1);
            aX_1(tmpIx) = NaN;
            aX_end(tmpIx) = NaN;
            
            
            aY_1 = (this.yPlanes(1) - this.sourcePoint(:,2)) ./ this.rayVec(:,2);
            aY_end = (this.yPlanes(end) - this.sourcePoint(:,2)) ./ this.rayVec(:,2);
            
            tmpIx = this.targetPoint(:,2) == this.sourcePoint(:,2);
            aY_1(tmpIx) = NaN;
            aY_end(tmpIx) = NaN;
            
            aZ_1 = (this.zPlanes(1) - this.sourcePoint(:,3)) ./ this.rayVec(:,3);
            aZ_end = (this.zPlanes(end) - this.sourcePoint(:,3)) ./ this.rayVec(:,3);
            
            tmpIx = this.targetPoint(:,3) == this.sourcePoint(:,3);
            aZ_1(tmpIx) = NaN;
            aZ_end(tmpIx) = NaN;
            
            
            % eq 5
            % Compute the \alpha_{min} and \alpha_{max} in terms of parametric values
            % given by equation 4.
            alpha_limits(:,1) = max([zeros(nRays,1) min(aX_1,aX_end) min(aY_1,aY_end) min(aZ_1,aZ_end)],[],2);
            alpha_limits(:,2) = min([ones(nRays,1) max(aX_1,aX_end) max(aY_1,aY_end) max(aZ_1,aZ_end)],[],2);
            
            % eq 6
            % Calculate the range of indeces who gives parametric values for
            % intersected planes.
            
            % get indices
            ixTpBiggerSp = this.targetPoint > this.sourcePoint; %this.rayVec > 0
            ixTpSmallerSp = this.targetPoint < this.sourcePoint; %this.rayVec < 0
            
            
            
            alphaTmp = NaN(nRays,2);
            alphaTmp(ixTpBiggerSp(:,1),1) = alpha_limits(ixTpBiggerSp(:,1),1);
            alphaTmp(ixTpSmallerSp(:,1),1) = alpha_limits(ixTpSmallerSp(:,1),2);
            alphaTmp(ixTpBiggerSp(:,1),2) = alpha_limits(ixTpBiggerSp(:,1),2);
            alphaTmp(ixTpSmallerSp(:,1),2) = alpha_limits(ixTpSmallerSp(:,1),1);
            
            i_min = this.numPlanes(1) - (this.xPlanes(end) - alphaTmp(:,1) .* this.rayVec(:,1) - this.sourcePoint(:,1))./this.resolution.x;
            i_max = 1          + (this.sourcePoint(:,1) + alphaTmp(:,2) .* this.rayVec(:,1) - this.xPlanes(1))./this.resolution.x;
            %rounding
            i_min = ceil(1/1000 * (round(1000*i_min)));
            i_max = floor(1/1000 * (round(1000*i_max)));
            
            
            %j
            alphaTmp(:) = NaN;
            alphaTmp(ixTpBiggerSp(:,2),1) = alpha_limits(ixTpBiggerSp(:,2),1);
            alphaTmp(ixTpSmallerSp(:,2),1) = alpha_limits(ixTpSmallerSp(:,2),2);
            alphaTmp(ixTpBiggerSp(:,2),2) = alpha_limits(ixTpBiggerSp(:,2),2);
            alphaTmp(ixTpSmallerSp(:,2),2) = alpha_limits(ixTpSmallerSp(:,2),1);
            
            j_min = this.numPlanes(2) - (this.yPlanes(end) - alphaTmp(:,1) .* this.rayVec(:,2) - this.sourcePoint(:,2))./this.resolution.y;
            j_max = 1          + (this.sourcePoint(:,2) + alphaTmp(:,2) .* this.rayVec(:,2) - this.yPlanes(1))./this.resolution.y;
            %rounding
            j_min = ceil(1/1000 * (round(1000*j_min)));
            j_max = floor(1/1000 * (round(1000*j_max)));
            
            %k
            alphaTmp(:) = NaN;
            alphaTmp(ixTpBiggerSp(:,3),1) = alpha_limits(ixTpBiggerSp(:,3),1);
            alphaTmp(ixTpSmallerSp(:,3),1) = alpha_limits(ixTpSmallerSp(:,3),2);
            alphaTmp(ixTpBiggerSp(:,3),2) = alpha_limits(ixTpBiggerSp(:,3),2);
            alphaTmp(ixTpSmallerSp(:,3),2) = alpha_limits(ixTpSmallerSp(:,3),1);
            
            k_min = this.numPlanes(3) - (this.zPlanes(end) - alphaTmp(:,1) .* this.rayVec(:,3) - this.sourcePoint(:,3))./this.resolution.z;
            k_max = 1          + (this.sourcePoint(:,3) + alphaTmp(:,2) .* this.rayVec(:,3) - this.zPlanes(1))./this.resolution.z;
            %rounding
            k_min = ceil(1/1000 * (round(1000*k_min)));
            k_max = floor(1/1000 * (round(1000*k_max)));
            
            % eq 7
            % For the given range of indices, calculate the paremetrics values who
            % represents intersections of the ray with the plane.
            
            planeIx = 1:length(this.xPlanes);
            planeGrid = repmat(this.xPlanes,nRays,1);
            planeGrid(planeIx < i_min | planeIx > i_max | (planeIx == i_min & planeIx == i_max) | isnan(i_min) | isnan(i_max)) = NaN;
            alpha_x = (planeGrid-this.sourcePoint(:,1))./(this.rayVec(:,1));
            
            planeIx = 1:length(this.yPlanes);
            planeGrid = repmat(this.yPlanes,nRays,1);
            planeGrid(planeIx < j_min | planeIx > j_max | (planeIx == j_min & planeIx == j_max) | isnan(j_min) | isnan(j_max)) = NaN;
            alpha_y = (planeGrid-this.sourcePoint(:,2))./(this.rayVec(:,2));
            
            planeIx = 1:length(this.zPlanes);
            planeGrid = repmat(this.zPlanes,nRays,1);
            planeGrid(planeIx < k_min | planeIx > k_max | (planeIx == k_min & planeIx == k_max) | isnan(k_min) | isnan(k_max)) = NaN;
            alpha_z = (planeGrid-this.sourcePoint(:,3))./(this.rayVec(:,3));
            
            % eq 8
            % Merge parametrics sets.
            % The following might look slow but is quite close to Matlab's
            % "unique" implementation
            alphas = sort([alpha_limits alpha_x alpha_y alpha_z],2); %NaN's are placed at the end when sorting in ascending order
            alphas(diff(alphas,1,2) == 0) = NaN; %Remove duplicates
            alphas = sort(alphas,2); %Again place NaN's at the end
            
            %Size Reduction (reduce NaN padding) for further computations
            maxNumColumns = max(sum(~isnan(alphas),2));
            alphas = alphas(:,1:maxNumColumns);
     end
       
        
        
        function this = initializeGeometry(this, ...
                resolution, ...
                cubeDim)
            %initializeGeometry Initializes CT geometry for Tracing
            %   Allows for vectorized input in source and target points
            % input
            %   resolution:     resolution of the cubes [mm/voxel]
            %   this.sourcePoint:    source point(s) of ray tracing. Either an 1x3
            %                   or Nx3 matrix (with N rays to trace). If
            %                   1x3 with multiple target points, the one
            %                   this.sourcePoint will be used for all rays.
            %   this.targetPoint:    target point(s) of ray tracing. Is an Nx3
            %                   matrix where N is the number of rays to
            %                   trace
            %   cubes:          cell array of cubes for ray tracing (it is possible to pass
            %                   multiple cubes for ray tracing to save computation time)
            %
            % output (see Siddon 1985 Medical Physics for a detailed description of the
            % variales)
            %   alphas          relative distance between start and endpoint for the
            %                    intersections with the cube
            %   l               lengths of intersestions with cubes
            %   rho             densities extracted from cubes
            %   d12             distance between start and endpoint of ray tracing
            %   ix              indices of hit voxels
            
            % Set private member variables for later use
            this.cubeDim = cubeDim;
            this.resolution = resolution;
            
            % Save the numbers of planes.            
            xNumPlanes = cubeDim(2) + 1;
            yNumPlanes = cubeDim(1) + 1;
            zNumPlanes = cubeDim(3) + 1;
            this.numPlanes = [xNumPlanes,yNumPlanes,zNumPlanes];
            
            % eq 3
            % Position of planes in millimeter. 0.5 because the central position
            % of the first voxel is at [resolution.x resolution.y resolution.z]
            
            this.xPlanes = resolution.x*(1:xNumPlanes) - .5*resolution.x;
            this.yPlanes = resolution.y*(1:yNumPlanes) - .5*resolution.y;
            this.zPlanes = resolution.z*(1:zNumPlanes) - .5*resolution.z;
            
            
            % Construct ray Vectors
            this.rayVec = this.targetPoint - this.sourcePoint;
            nRays = size(this.rayVec,1);
            
            
            
            
            
        end
        
        
        
    end
end