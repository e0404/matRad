function shapeInfo = tk_getSequencingParameters(seqResult,pln,stf,visBool)

mode = 'physical'; %'index' 'physical'
% physical: 81 leafes, 0mm = middle of leaf 41
bixelWidth = pln.bixelWidth;

% initializing variables
x_l = [];
x_r = [];
x_min = [];
x_max = [];
shapeIx = [];
shapeCounter = 0;
numOfShapes = ones(pln.numOfBeams,1);
leafIx = [];
lim_l = [];
lim_r = [];
IxOffset = 0;


if strcmp(mode,'index')
% loop over all beams
    for i=1:pln.numOfBeams

        % get number of shapes for each beam
        numOfShapes(i) = seqResult.beam(i).numOfShapes;

        % get x- and z-coordinates of bixels
        X = ones(stf(i).numOfRays,1)*NaN;
        Z = ones(stf(i).numOfRays,1)*NaN;

        for k=1:stf(i).numOfRays
            X(k) = stf(i).ray(k).rayPos_bev(:,1);
            Z(k) = stf(i).ray(k).rayPos_bev(:,3);
        end

        % create ray-map
        maxX = max(X);
        minX = min(X);

        maxZ = max(Z);
        minZ = min(Z);    

        dimOfRayMxX = (maxX-minX)/stf(i).bixelWidth + 1;
        dimOfRayMxZ = (maxZ-minZ)/stf(i).bixelWidth + 1;

        rayMx = zeros(dimOfRayMxZ,dimOfRayMxX);


        % get indices for x and z positions
        xPos = (X-minX)/stf(i).bixelWidth +1;
        zPos = (Z-minZ)/stf(i).bixelWidth +1;

        % get indices in the ray map
        IxInRayMX = zPos + (xPos-1)*dimOfRayMxZ;

        % fill ray-map
        rayMx(IxInRayMX) = 1;

        % fill map of bixel indices
        IxMaps(i).bixelIndices = NaN * ones(dimOfRayMxZ,dimOfRayMxX);
        counter = 1;

        for n=1:numel(IxMaps(i).bixelIndices)
            if ismember(n,IxInRayMX)
                IxMaps(i).bixelIndices(n) = counter;
                counter = counter + 1;
            end
        end
        IxMaps(i).bixelIndices = IxMaps(i).bixelIndices + IxOffset;
        IxOffset = IxOffset + numel(X);

        % get leaf number
        leafNum = Z/bixelWidth + 41;
        leafNum = unique(leafNum);% only keep each leaf once

        % loop over all shapes
        for j=1:seqResult.beam(i).numOfShapes
            shapeCounter = shapeCounter +1; % get index of current shape

            if j==1
                % get the limits for the leaf pairs
                [shapeLim_l, shapeLim_r] = tk_getLeafPos(rayMx);
            end


            tempShape = seqResult.beam(i).shapes(:,:,j);
            [Ix_l, Ix_r, minIx, maxIx] = tk_getLeafPos(tempShape);

            %visualize shape
            if visBool
                figure
                titleString = ['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(j)];
                tk_visSingleShape(tempShape,Ix_l,Ix_r,titleString)
                title(sprintf(['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(j)],'Fontsize',14))
                close all
            end

            % add leaf positions to position vectors
            x_l = [x_l; Ix_l];
            x_r = [x_r; Ix_r];
            x_min = [x_min; minIx];
            x_max = [x_max; maxIx];
            shapeIx = [shapeIx; shapeCounter*ones(size(Ix_l,1),1)];
            leafIx = [leafIx; leafNum];
            lim_l = [lim_l; shapeLim_l];
            lim_r = [lim_r; shapeLim_r];      
        end

    end

    shapeInfo.x_l = x_l;
    shapeInfo.x_r = x_r;
    shapeInfo.x_min = x_min;
    shapeInfo.x_max = x_max;
    shapeInfo.leafIx = leafIx;
    shapeInfo.shapeIx = shapeIx;
    shapeInfo.numOfShapes = numOfShapes;
    shapeInfo.lim_l = lim_l;
    shapeInfo.lim_r = lim_r;
    shapeInfo.IxMaps = IxMaps;
end

if strcmp(mode,'physical') % NOCH NICHT GEMACHT!!!!NUR KOPIE!!!
    % loop over all beams
    for i=1:pln.numOfBeams
        % get number of shapes for each beam
        numOfShapes(i) = seqResult.beam(i).numOfShapes;

        % get x- and z-coordinates of bixels
        X = ones(stf(i).numOfRays,1)*NaN;
        Z = ones(stf(i).numOfRays,1)*NaN;

        for k=1:stf(i).numOfRays
            X(k) = stf(i).ray(k).rayPos_bev(:,1);
            Z(k) = stf(i).ray(k).rayPos_bev(:,3);
        end

        % get leaf number
        leafNum = Z/bixelWidth + 41;
        leafNum = unique(leafNum);% only keep each leaf once
        
                % create ray-map
                maxX = max(X);
                minX = min(X);

                maxZ = max(Z);
                minZ = min(Z);    

                dimOfRayMxX = (maxX-minX)/stf(i).bixelWidth + 1;
                dimOfRayMxZ = (maxZ-minZ)/stf(i).bixelWidth + 1;

                rayMx = zeros(dimOfRayMxZ,dimOfRayMxX);

                % get indices for x and z positions
                xPos = (X-minX)/stf(i).bixelWidth +1;
                zPos = (Z-minZ)/stf(i).bixelWidth +1;

                % get indices in the ray map
                IxInRayMX = zPos + (xPos-1)*dimOfRayMxZ;

                % fill ray-map
                rayMx(IxInRayMX) = 1;

                % fill map of bixel indices
                IxMaps(i).bixelIndices = NaN * ones(dimOfRayMxZ,dimOfRayMxX);
                counter = 1;

                for n=1:numel(IxMaps(i).bixelIndices)
                    if ismember(n,IxInRayMX)
                        IxMaps(i).bixelIndices(n) = counter;
                        counter = counter + 1;
                    end
                end
                IxMaps(i).bixelIndices = IxMaps(i).bixelIndices + IxOffset;
                IxOffset = IxOffset + numel(X);
        
        for j=1:seqResult.beam(i).numOfShapes
            shapeCounter = shapeCounter +1; % get index of current shape

            if j==1
                % get the limits for the leaf pairs
                [shapeLim_l, shapeLim_r] = tk_getLeafPos(rayMx);
            end


            tempShape = seqResult.beam(i).shapes(:,:,j);
            [Ix_l, Ix_r, minIx, maxIx] = tk_getLeafPos(tempShape);

            lineIx = [];
            for o=1:numel(Ix_l)
                lineIx = [lineIx; o];
            end

bixelIx = IxMaps(i).bixelIndices(lineIx,Ix_l+1);
            
            %visualize shape
            if visBool
                figure
                titleString = ['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(j)];
                tk_visSingleShape(tempShape,Ix_l,Ix_r,titleString)
                title(sprintf(['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(j)],'Fontsize',14))
                close all
            end

            
            % calculate physical leaf positions
                % 1. get bixel indices from leaf index
                
            
            % add leaf positions to position vectors
            x_l = [x_l; Ix_l];
            x_r = [x_r; Ix_r];
            x_min = [x_min; minIx];
            x_max = [x_max; maxIx];
            shapeIx = [shapeIx; shapeCounter*ones(size(Ix_l,1),1)];
            leafIx = [leafIx; leafNum];
%             lim_l = [lim_l; shapeLim_l];
%             lim_r = [lim_r; shapeLim_r];      
        end
    end
    
    shapeInfo.x_l = x_l;
    shapeInfo.x_r = x_r;
    shapeInfo.x_min = x_min;
    shapeInfo.x_max = x_max;
    shapeInfo.leafIx = leafIx;
    shapeInfo.shapeIx = shapeIx;
    shapeInfo.numOfShapes = numOfShapes;
%     shapeInfo.lim_l = lim_l;
%     shapeInfo.lim_r = lim_r;
%     shapeInfo.IxMaps = IxMaps;
    
end

    
end