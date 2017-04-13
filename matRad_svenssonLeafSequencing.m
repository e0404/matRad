function resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,numToKeep,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multileaf collimator leaf sequencing algorithm for intensity modulated
% beams with multiple static segments according to Siochi (1999)
% International Journal of Radiation Oncology * Biology * Physics,
% originally implemented in PLUNC (https://sites.google.com/site/planunc/)
%
% Implented in matRad by Eric Christiansen, Emily Heath, and Tong Xu
%
% call
%   resultSequencing =
%   matRad_siochiLeafSequencing(w,stf,pln,numOfLevels,visBool)
%
% input
%   resultGUI:          resultGUI struct to which the output data will be
%   added, if
%                       this field is empty resultGUI struct will be
%                       created
%   stf:                matRad steering information struct
%   numOfLevels:        number of stratification levels
%   visBool:            toggle on/off visualization (optional)
%
% output
%   resultGUI:          matRad result struct containing the new dose cube
%   as well as the corresponding weights
%
% References
%   [1] https://www.ncbi.nlm.nih.gov/pubmed/10078655
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if visBool not set toogle off visualization
if nargin < 5
    numToKeep = 0;
elseif nargin < 6
     visBool = 0;
end

if pln.VMAT
    %First beam sweeps right-to-left, next left-to-right, ...
    inversion = 1;
else
    if numToKeep == 0
        error('\n\nIf not doing VMAT, must keep a non-zero number of apertures.\n\n');
    end
end
    

numOfBeams = numel(stf);

if visBool
    % create the sequencing figure
    sz = [800 1000]; % figure size
    screensize = get(0,'ScreenSize');
    xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
    ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically
    seqFig = figure('position',[xpos,ypos,sz(2),sz(1)]);     
end

offset = 0;


for i = 1:numOfBeams
    numOfRaysPerBeam = stf(i).numOfRays;
    if isfield(stf(i),'initializeBeam') && ~stf(i).initializeBeam
        sequencing.w(1+offset:numOfRaysPerBeam+offset,1) = 0;
        offset = offset + numOfRaysPerBeam;
        if ~pln.VMAT
            sequencing.beam(i).numOfShapes = 0;
        end
        continue %if this is not a beam to be initialized, continue to next iteration without generating segments
    end
    
    % get relevant weights for current beam
    wOfCurrBeams = resultGUI.wUnsequenced(1+offset:numOfRaysPerBeam+offset);%REVIEW OFFSET
    
    X = ones(numOfRaysPerBeam,1)*NaN;
    Z = ones(numOfRaysPerBeam,1)*NaN;
    
    for segment=1:stf(i).numOfRays
        X(segment) = stf(i).ray(segment).rayPos_bev(:,1);
        Z(segment) = stf(i).ray(segment).rayPos_bev(:,3);
    end
    
    % sort bixels into matrix
    minX = min(X);
    maxX = max(X);
    minZ = min(Z);
    maxZ = max(Z);
    
    dimOfFluenceMxX = (maxX-minX)/stf(i).bixelWidth + 1;
    dimOfFluenceMxZ = (maxZ-minZ)/stf(i).bixelWidth + 1;
    
    %Create the fluence matrix.
    fluenceMx = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    
    % Calculate X and Z position of every fluence's matrix spot z axis =
    % axis of leaf movement!
    xPos = (X-minX)/stf(i).bixelWidth+1;
    zPos = (Z-minZ)/stf(i).bixelWidth+1;
    
    % Make subscripts for fluence matrix
    indInFluenceMx = zPos + (xPos-1)*dimOfFluenceMxZ;
    
    %Save weights in fluence matrix.
    fluenceMx(indInFluenceMx) = wOfCurrBeams;
    
    t = fluenceMx/(pln.defaultDoseRate/dij.weightToMU);
    [dtdx,~] = gradient(t);
    dtdx = dtdx./stf(i).bixelWidth;
    
    dtdx_pos = dtdx;
    dtdx_neg = dtdx;
    dtdx_pos(dtdx_pos < 0) = 0;
    dtdx_neg(dtdx_neg >= 0) = 0;
    
    tL = zeros(size(t,1),size(t,2)+1);
    tR = zeros(size(t,1),size(t,2)+1);
    
    tMin = zeros(size(t,1),1);
    tGrad = zeros(size(t,1),1);
    
    %first, find longest trajectory time
    for row = 1:size(tL,1)
        tMin(row) = (size(tL,2)-1)*stf(i).bixelWidth/pln.defaultLeafSpeed+stf(i).bixelWidth*sum(dtdx_pos(row,1:(size(tL,2)-1)));
        tGrad(row) = stf(i).bixelWidth*sum(dtdx_pos(row,1:(size(tL,2)-1)));
    end
    tMinMax = max(tMin);
    
    maxLeafSpeed = (size(tL,2)-1)*stf(i).bixelWidth./(tMinMax-tGrad);
    
    for row = 1:size(tL,1)
        for col = 1:size(tL,2)
            tL(row,col) = (col-1)*stf(i).bixelWidth/maxLeafSpeed(row)+stf(i).bixelWidth*sum(dtdx_pos(row,1:(col-1)));
            tR(row,col) = (col-1)*stf(i).bixelWidth/maxLeafSpeed(row)-stf(i).bixelWidth*sum(dtdx_neg(row,1:(col-1)));
        end
    end
    
    %{
    for row = 1:size(t,1)
        tL(row,:) = tL(row,:)*tMinMax./tMin(row);
        tR(row,:) = tR(row,:)*tMinMax./tMin(row);
    end
    %}
    
    if numToKeep ~= 0 && ~pln.VMAT
        k = numToKeep;
        shapes = zeros(dimOfFluenceMxZ,dimOfFluenceMxX,k);
        shapesWeight = zeros(k,1);
        
        %tStep = tMinMax/(k-1);
        %tSample = 0;
        
        tStep = tMinMax*ones(k,1)/k;
        tSample = cumsum(tStep)-tStep/2;
        for shape = 1:k
            for row = 1:size(tL,1)
                xL = round(interp1(tL(row,:),(1:size(tL,2)),tSample(shape),'linear','extrap'));
                xR = round(interp1(tR(row,:),(1:size(tL,2))-1,tSample(shape),'linear','extrap'));
                
                xL = min(xL,size(tL,2));
                xR = min(xR,size(tL,2)-1);
                
                shapes(row,xL:xR,shape) = 1;
            end
            shapesWeight(shape) = (pln.defaultDoseRate/dij.weightToMU)*tStep(shape);
        end
        
        sequencing.beam(i).numOfShapes  = k;
        sequencing.beam(i).shapes       = shapes(:,:,1:k);
        sequencing.beam(i).shapesWeight = shapesWeight(1:k);
        sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
        %sequencing.beam(i).fluence      = D_0;
        sequencing.beam(i).sum          = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
        
        for segment = 1:sequencing.beam(i).numOfShapes
            sequencing.beam(i).sum = sequencing.beam(i).sum+sequencing.beam(i).shapes(:,:,segment)*sequencing.beam(i).shapesWeight(segment);
        end
        sequencing.w(1+offset:numOfRaysPerBeam+offset,1) = sequencing.beam(i).sum(indInFluenceMx);
        
        offset = offset + numOfRaysPerBeam;
        
    elseif pln.VMAT
        %Divide segments into n sectors, sample leaf trajectory at n points
        %to get left and right leaf positions.  Spread apertures to
        %children of the initialized beam angle.  Make assumption that
        %gantry rotation speed and dose rate are constant over each arc
        %sector.
        
        sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
        
        numSectors = stf(i).numOfBeamChildren;
        childrenAngles = stf(i).beamChildrenGantryAngles;
        childrenIndex = stf(i).beamChildrenIndex;
        borderAngles = stf(i).borderAngles;
        
        angleStep = diff(childrenAngles);
        angleStep = padarray(angleStep,[1,0],angleStep(end),'pre');
        tStep = tMinMax*angleStep/sum(angleStep);
        tSample = cumsum(tStep)-tStep/2;
        
        sector = 1;
        inversion = -1*inversion; % =1 (-1) for even (odd) init gantry angles
        sector0 = (inversion*(1-numSectors)+(1+numSectors))/2; % =1 (k) for even (odd) init gantry angles
        sector1 = numSectors+1-sector0;% =k (1) for even (odd) init gantry angles
        
        for segment = sector0:inversion:sector1 % = 1-to-k in forward (reverse) order for even (odd) values of i
            shape = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
            
            for row = 1:size(tL,1)
                xL = round(interp1(tL(row,:),(1:size(tL,2)),tSample(segment),'linear','extrap'));
                xR = round(interp1(tR(row,:),(1:size(tL,2))-1,tSample(segment),'linear','extrap'));
                
                xL = min(xL,size(tL,2));
                xR = min(xR,size(tL,2)-1);
                
                shape(row,xL:xR) = 1;
            end
            shapeWeight = (pln.defaultDoseRate/dij.weightToMU)*tStep(segment);
            
            childIndex = childrenIndex(sector);
            
            sequencing.beam(childIndex).numOfShapes = 1;
            sequencing.beam(childIndex).shapes = shape;
            sequencing.beam(childIndex).shapesWeight = shapeWeight;
            
            big = max([i-1 childIndex-1]);
            small = min([i childIndex]);
            numOfRaysBN = sum([stf(small:big).numOfRays]); %num of rays between current beam and child
            inv = (childIndex-i)./(big-small+1); % = 1 (-1) if childIndex is bigger (smaller) than i
            sequencing.beam(childIndex).bixelIx = sequencing.beam(i).bixelIx+inv*numOfRaysBN;
            
            sequencing.beam(childIndex).fluence = shape;
            sequencing.beam(childIndex).sum = shapeWeight*shape;
            
            sector = sector+1;
        end
        
        offset = offset + numOfRaysPerBeam;
    end
end

if pln.VMAT
    gantryAngles = [stf.gantryAngle];
    optGantryAngles = [stf([stf.optimizeBeam]).gantryAngle];
    
    for i = 1:numel(optGantryAngles)
        currInd = find(gantryAngles==optGantryAngles(i));
        
        if stf(currInd).initializeBeam
            
            for j = 1:stf(currInd).numOfBeamSubChildren
                %Prevents matRad_sequencing2ApertureInfo from attempting to
                %convert shape to aperturevec for subchildren
                sequencing.beam(stf(currInd).beamSubChildrenIndex(j)).numOfShapes = 0;
            end
        end
        
        if i ~= numel(optGantryAngles)
            sequencing.beam(currInd).gantryRot = pln.defaultGantryRot; %gantry rotation rate until next opt angle
            nextInd = find(gantryAngles==optGantryAngles(i+1));
            sequencing.beam(currInd).MURate = dij.weightToMU.*sequencing.beam(currInd).shapesWeight.*sequencing.beam(currInd).gantryRot./(stf(nextInd).gantryAngle-stf(currInd).gantryAngle); %dose rate until next opt angle
            %Rescale weight to represent only this beam; assume piecewise
            %linear doserate
            sequencing.beam(currInd).shapesWeight = sequencing.beam(currInd).shapesWeight.*((stf(currInd+1).gantryAngle-stf(currInd).gantryAngle)./(stf(nextInd).gantryAngle-stf(currInd).gantryAngle));
        else
            sequencing.beam(currInd).gantryRot = pln.defaultGantryRot; %gantry rotation rate until next opt angle
            prevInd = find(gantryAngles==optGantryAngles(i-1));
            sequencing.beam(currInd).MURate = dij.weightToMU.*sequencing.beam(currInd).shapesWeight.*sequencing.beam(currInd).gantryRot./(stf(currInd).gantryAngle-stf(prevInd).gantryAngle); %dose rate until next opt angle
            %Rescale weight to represent only this beam; assume piecewise
            %linear doserate
            if currInd+1 <= numel(gantryAngles)
                sequencing.beam(currInd).shapesWeight = sequencing.beam(currInd).shapesWeight.*((stf(currInd+1).gantryAngle-stf(currInd).gantryAngle)./(stf(currInd).gantryAngle-stf(prevInd).gantryAngle));
            else
                sequencing.beam(currInd).shapesWeight = sequencing.beam(currInd).shapesWeight.*((stf(currInd).gantryAngle-stf(currInd-1).gantryAngle)./(stf(currInd).gantryAngle-stf(prevInd).gantryAngle));
            end
        end
        
    end
    
    %Calculate w using matRad functions
    sequencing.weightToMU = dij.weightToMU;
    resultGUI.apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf,pln.VMAT);
    
    
    %matRad_daoVec2ApertureInfo will interpolate subchildren gantry
    %segments
    resultGUI.apertureInfo = matRad_daoVec2ApertureInfo(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector,1);
    
    % LEAF TRAVEL / DEGREE
    
    %calculate max leaf speed
    resultGUI.apertureInfo = matRad_maxLeafSpeed(resultGUI.apertureInfo);
    
    %optimize delivery
    resultGUI = matRad_optDelivery(resultGUI,pln,0);
    
    sequencing.w = resultGUI.apertureInfo.bixelWeights;
    
else
    sequencing.weightToMU = dij.weightToMU;
    resultGUI.apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf);
end

resultGUI.w          = sequencing.w;
resultGUI.wSequenced = sequencing.w;

resultGUI.sequencing   = sequencing;

%EXAMINE
resultGUI.physicalDose = reshape(dij.physicalDose{1} * sequencing.w,dij.dimensions);

% if weights exists from an former DAO remove it
if isfield(resultGUI,'wDao')
    resultGUI = rmfield(resultGUI,'wDao');
end





end
