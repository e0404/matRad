function shapeInfo = tk_getParameters(Sequencing,stf,pln,visBool)

% MLC parameters:
bixelWidth = pln.bixelWidth; % [mm]
numOfMLCLeafPairs = 80;
%     define central leaf pair (here we want the 0mm position to be in the
%     center of a leaf pair (e.g. leaf 41 stretches from -2.5mm to 2.5mm
%     for a bixel/leafWidth of 5mm)
centralLeafPair = floor(median(1:numOfMLCLeafPairs));

% initializing variables
bixelIndOffset = 0; % used for creation of bixel index maps
totalNumOfBixels = sum([stf(:).totalNumOfBixels]);
vectorIndex = totalNumOfBixels + 1; % used for bookkeeping in the vector for optimization





% loop over all beams
for i=1:pln.numOfBeams
    
    %% 1. read stf and create maps (Ray & Bixelindex)
    
    % get x- and z-coordinates of bixels
    X = ones(stf(i).numOfRays,1)*NaN;
    Z = ones(stf(i).numOfRays,1)*NaN;
    for j=1:stf(i).numOfRays
        X(j) = stf(i).ray(j).rayPos_bev(:,1);
        Z(j) = stf(i).ray(j).rayPos_bev(:,3);
    end
    
    % create ray-map
        maxX = max(X);
        minX = min(X);    
        maxZ = max(Z);
        minZ = min(Z);

        dimX = (maxX-minX)/stf(i).bixelWidth + 1;
        dimZ = (maxZ-minZ)/stf(i).bixelWidth + 1;

        rayMap = zeros(dimZ,dimX);

        % get indices for x and z positions
        xPos = (X-minX)/stf(i).bixelWidth +1;
        zPos = (Z-minZ)/stf(i).bixelWidth +1;

        % get indices in the ray-map
        indInRay = zPos + (xPos-1)*dimZ;

        % fill ray-map
        rayMap(indInRay) = 1;
    
    % create map of bixel indices
        bixelIndMap = NaN * ones(dimZ,dimX);
        counter = 1;

        for k=1:numel(bixelIndMap)
            if ismember(k,indInRay)
                bixelIndMap(k) = counter;
                counter = counter + 1;
            end
        end
        bixelIndMap = bixelIndMap + bixelIndOffset;
        bixelIndOffset = bixelIndOffset + numel(X);
    % store physical position of first entry in bixelIndMap
        posOfCornerBixel = [minX 0 minZ];
    
    %% 2. Get leaf limits
    % leaf limits can be extracted from the leaf map
        % initializing limits
        lim_l = NaN * ones(dimZ,1);
        lim_r = NaN * ones(dimZ,1);
        % looping oder leaf pairs
        for l=1:dimZ
            lim_lInd = find(rayMap(l,:),1,'first');
            lim_rInd = find(rayMap(l,:),1,'last');
            % the physical position [mm] can be calculated from the indices
            lim_l(l) = (lim_lInd-1)*bixelWidth + minX - 1/2*bixelWidth;
            lim_r(l) = (lim_rInd-1)*bixelWidth + minX + 1/2*bixelWidth;
        end
        
    %% 3. Get leaf positions
    % leaf positions can be extracted from the shapes created in Sequencing
        for m=1:Sequencing.beam(i).numOfShapes
            % loading shape from Sequencing result
            shapeMap = Sequencing.beam(i).shapes(:,:,m);
            % get left and right leaf index from shapemap
                % initializing limits
                leftLeafPos = NaN * ones(dimZ,1);
                rightLeafPos = NaN * ones(dimZ,1);
                % looping over leaf pairs
                for l=1:dimZ
                    try %try to find open bixel
                        leftLeafPosInd = find(shapeMap(l,:),1,'first');
                        rightLeafPosInd = find(shapeMap(l,:),1,'last');
                        % the physical position [mm] can be calculated from the indices
                        leftLeafPos(l) = (leftLeafPosInd-1)*bixelWidth...
                                        + minX - 1/2*bixelWidth;
                        rightLeafPos(l) = (rightLeafPosInd-1)*bixelWidth...
                                        + minX + 1/2*bixelWidth;
                    catch % if no bixel is open, use limits from Ray positions
                        leftLeafPos(l) = (lim_l(l)+lim_r(l))/2;
                        rightLeafPos(l) = leftLeafPos(l);
                    end
                end
            % save data for each shape of this beam
            shapeInfo.beam(i).shape(m).leftLeafPos = leftLeafPos;
            shapeInfo.beam(i).shape(m).rightLeafPos = rightLeafPos;
            shapeInfo.beam(i).shape(m).weight = Sequencing.beam(i).shapesWeight(m);
            shapeInfo.beam(i).shape(m).shapeMap = shapeMap;
            shapeInfo.beam(i).shape(m).index = vectorIndex;
            
            %update index for bookkeeping
            vectorIndex = vectorIndex + dimZ;
            
            %visualize result
            if visBool
                figure
                titleString = ['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(m)];
                tk_visSingleShape(shapeMap,leftLeafPos,rightLeafPos,...
                    titleString,minX-bixelWidth/2,maxX+bixelWidth/2)
                title(sprintf(['BeamNr: ' int2str(i) ' ShapeNr: ' int2str(m)],'Fontsize',14))
                
                pause(0.5)
                close all
            end
        end
        
    %% 4. Get z-coordinates of active leaf pairs        
%     get z-coordinates from bixel positions
        leafPairPos = unique(Z);
        
%     find upmost and downmost leaf pair
        topLeafPairPos = maxZ;
        bottomLeafPairPos = minZ;
        
        topLeafPair = centralLeafPair - topLeafPairPos/bixelWidth;
        bottomLeafPair = centralLeafPair - bottomLeafPairPos/bixelWidth;
        
%     create bool map of active leaf pairs
        isActiveLeafPair = zeros(numOfMLCLeafPairs,1);
        isActiveLeafPair(topLeafPair:bottomLeafPair) = 1;
        
    %% 5. create MLC window
%       getting the dimensions of the MLC in order to be able to plot the
%       shapes using physical coordinates
        MLCWindow = [minX-bixelWidth/2 maxX+bixelWidth/2 ...
                    minZ-bixelWidth/2 maxZ+bixelWidth/2];
           
    
    
    
    %% save data for each beam
    shapeInfo.beam(i).numOfShapes = Sequencing.beam(i).numOfShapes;
    shapeInfo.beam(i).numOfActiveLeafPairs = dimZ;
    shapeInfo.beam(i).leafPairPos = leafPairPos;
    shapeInfo.beam(i).isActiveLeafPair = isActiveLeafPair;
    shapeInfo.beam(i).centralLeafPair = centralLeafPair;
    shapeInfo.beam(i).lim_l = lim_l;
    shapeInfo.beam(i).lim_r = lim_r;
    shapeInfo.beam(i).bixelIndMap = bixelIndMap;
    shapeInfo.beam(i).posOfCornerBixel = posOfCornerBixel;
    shapeInfo.beam(i).MLCWindow = MLCWindow;
end




%% save global data
shapeInfo.bixelWidth = bixelWidth;
shapeInfo.numOfMLCLeafPairs = numOfMLCLeafPairs;



end