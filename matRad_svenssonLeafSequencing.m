function resultGUI = matRad_svenssonLeafSequencing(resultGUI,stf,dij,pln,visBool)
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
    visBool = 0;
end

if pln.VMAT
    %First beam sweeps right-to-left, next left-to-right, ...
    leafDir = 1;
end
sequencing.VMAT = pln.VMAT;
sequencing.dynamic = pln.dynamic;

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
    
    temp = zeros(size(fluenceMx));
    for row = 1:dimOfFluenceMxZ
        temp(row,:) = imgaussfilt(fluenceMx(row,:),1);
    end
    fluenceMx = temp;
    clear temp
    
    t = fluenceMx/(pln.defaultDoseRate/dij.weightToMU);
    
    minInd = zeros(dimOfFluenceMxZ,1);
    maxInd = zeros(dimOfFluenceMxZ,1);
    
    for row = 1:dimOfFluenceMxZ
        minInd(row) = find(t(row,:) > 0,1,'first');
        maxInd(row) = find(t(row,:) > 0,1,'last');
    end
    
    tpad = padarray(t,[0 1],0,'both');
    %dtdx = grad(t)./stf(i).bixelWidth;
    
    dtdx = diff(tpad,1,2)./stf(i).bixelWidth;
    dtdx_pos = dtdx;
    dtdx_neg = dtdx;
    dtdx_pos(dtdx_pos < 0) = 0;
    dtdx_neg(dtdx_neg >= 0) = 0;
    
    tL = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    tR = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    
    %first, find longest trajectory time
    tMin = (maxInd-minInd+1)*stf(i).bixelWidth/pln.defaultLeafSpeed+stf(i).bixelWidth*sum(dtdx_pos,2);
    tGrad = stf(i).bixelWidth*sum(dtdx_pos,2);
    %{
    for row = 1:size(tL,1)
        tMin(row) = (size(tL,2)-1)*stf(i).bixelWidth/pln.defaultLeafSpeed+stf(i).bixelWidth*sum(dtdx_pos(row,1:(size(tL,2)-1)));
        tGrad(row) = stf(i).bixelWidth*sum(dtdx_pos(row,1:(size(tL,2)-1)));
    end
    %}
    tMinMax = max(tMin);
    
    maxLeafSpeed = (maxInd-minInd+1)*stf(i).bixelWidth./(tMinMax-tGrad);
    
    for row = 1:dimOfFluenceMxZ
        for col = 1:dimOfFluenceMxX
            tL(row,col) = (col-minInd(row)+1)*stf(i).bixelWidth/maxLeafSpeed(row)+stf(i).bixelWidth*sum(dtdx_pos(row,1:col));
            tR(row,col) = (col-minInd(row)+1)*stf(i).bixelWidth/maxLeafSpeed(row)-stf(i).bixelWidth*sum(dtdx_neg(row,1:col));
        end
    end
    
    %{
    for row = 1:size(t,1)
        tL(row,:) = tL(row,:)*tMinMax./tMin(row);
        tR(row,:) = tR(row,:)*tMinMax./tMin(row);
    end
    %}
    
    if ~pln.VMAT
        k = pln.numApertures;
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
        
    elseif pln.VMAT && pln.dynamic
        %Divide segments into n sectors, sample leaf trajectory at n points
        %to get left and right leaf positions.  Spread apertures to
        %children of the initialized beam angle.  Make assumption that
        %gantry rotation speed and dose rate are constant over each arc
        %sector.
        
        sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
        
        angleInit_I = stf(i).initAngleBorders(1);
        angleInit_F = stf(i).initAngleBorders(2);
        totAngleInit = stf(i).initAngleBordersDiff;
        numSectors = stf(i).numOfBeamChildren;
        childrenIndex = stf(i).beamChildrenIndex;
        
        leafDir = -1*leafDir; % =1 (-1) for even (odd) init gantry angles
        xPos = unique(X');
        if leafDir == -1
            tL = tMinMax-tL;
            tR = tMinMax-tR;
        end
        
        for segment = 1:numSectors
            childIndex = childrenIndex(segment);
            
            angleSample_I = stf(childIndex).doseAngleBorders(1);
            angleSample_F = stf(childIndex).doseAngleBorders(2);
            angleSample = stf(childIndex).gantryAngle;
            
            tSample_I = (angleSample_I-angleInit_I)*tMinMax./totAngleInit;
            tSample_F = (angleSample_F-angleInit_I)*tMinMax./totAngleInit;
            tSample = (angleSample-angleInit_I)*tMinMax./totAngleInit;
            
            xL_I = zeros(size(tL,1),1);
            xR_I = zeros(size(tL,1),1);
            xL_F = zeros(size(tL,1),1);
            xR_F = zeros(size(tL,1),1);
            xL = zeros(size(tL,1),1);
            xR = zeros(size(tL,1),1);
            
            for row = 1:size(tL,1)
                xL_I(row) = interp1(tL(row,:),xPos,tSample_I,'linear','extrap');
                xR_I(row) = interp1(tR(row,:),xPos,tSample_I,'linear','extrap');
                
                xL_F(row) = interp1(tL(row,:),xPos,tSample_F,'linear','extrap');
                xR_F(row) = interp1(tR(row,:),xPos,tSample_F,'linear','extrap');
                
                xL(row) = interp1(tL(row,:),xPos,tSample,'linear','extrap');
                xR(row) = interp1(tR(row,:),xPos,tSample,'linear','extrap');
            end
            shapeWeight = (pln.defaultDoseRate/dij.weightToMU)*(tSample_F-tSample_I);
            
            sequencing.beam(childIndex).shape(1).leftLeafPos_I = xL_I;
            sequencing.beam(childIndex).shape(1).rightLeafPos_I = xR_I;
            sequencing.beam(childIndex).shape(1).leftLeafPos_F = xL_F;
            sequencing.beam(childIndex).shape(1).rightLeafPos_F = xR_F;
            sequencing.beam(childIndex).shape(1).leftLeafPos = xL;
            sequencing.beam(childIndex).shape(1).rightLeafPos = xR;
            
            sequencing.beam(childIndex).numOfShapes = 1;
            sequencing.beam(childIndex).leafDir = leafDir;
            
            sequencing.beam(childIndex).shapesWeight = shapeWeight;
            
            big = max([i-1 childIndex-1]);
            small = min([i childIndex]);
            numOfRaysBN = sum([stf(small:big).numOfRays]); %num of rays between current beam and child
            inv = (childIndex-i)./(big-small+1); % = 1 (-1) if childIndex is bigger (smaller) than i
            sequencing.beam(childIndex).bixelIx = sequencing.beam(i).bixelIx+inv*numOfRaysBN;
        end
        
        offset = offset + numOfRaysPerBeam;
    end
end

if pln.VMAT
    gantryAngles = [stf.gantryAngle];
    optGantryAngles = [stf([stf.optimizeBeam]).gantryAngle];
    
    for i = 1:numel(optGantryAngles)
        optInd = find(gantryAngles==optGantryAngles(i));
        
        if stf(optInd).initializeBeam
            for j = 1:stf(optInd).numOfBeamSubChildren
                %Prevents matRad_sequencing2ApertureInfo from attempting to
                %convert shape to aperturevec for subchildren
                sequencing.beam(stf(optInd).beamSubChildrenIndex(j)).numOfShapes = 0;
            end
        end
        
        sequencing.beam(optInd).gantryRot = pln.defaultGantryRot; %gantry rotation rate until next opt angle
        sequencing.beam(optInd).MURate = dij.weightToMU.*sequencing.beam(optInd).shapesWeight.*sequencing.beam(optInd).gantryRot./diff(stf(optInd).optAngleBorders); %dose rate until next opt angle
        %Rescale weight to represent only this control point; weight will be shared
        %with the interpolared control points in matRad_daoVec2ApertureInfo
        sequencing.beam(optInd).shapesWeight = sequencing.beam(optInd).shapesWeight.*stf(optInd).timeFacCurr;
        
    end
    
    %Calculate w using matRad functions
    sequencing.weightToMU = dij.weightToMU;
    sequencing.jacobi = pln.jacobi;
    resultGUI.apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf);
    
    
    %matRad_daoVec2ApertureInfo will interpolate subchildren gantry
    %segments
    resultGUI.apertureInfo.updateJacobi = true;
    if pln.dynamic
        resultGUI.apertureInfo = matRad_daoVec2ApertureInfo_VMATdynamic(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);
    else
        resultGUI.apertureInfo = matRad_daoVec2ApertureInfo_VMATstatic(resultGUI.apertureInfo,resultGUI.apertureInfo.apertureVector);
    end
    resultGUI.apertureInfo.updateJacobi = false;
    
    % LEAF TRAVEL / DEGREE
    
    %calculate max leaf speed
    resultGUI.apertureInfo = matRad_maxLeafSpeed(resultGUI.apertureInfo);
    
    %optimize delivery
    resultGUI = matRad_optDelivery(resultGUI,pln,0);
    resultGUI.apertureInfo = matRad_maxLeafSpeed(resultGUI.apertureInfo);
    
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


function g = grad(f)
% f is a m x n matrix
% g is the gradient of f in the horizontal ([0 1]) direction

fpad = padarray(f,[0 1],0,'both');
g = zeros(size(fpad));

for i = 2:(size(fpad,2)-2)
    g(:,i) = (fpad(:,i+1)-fpad(:,i-1))/2;
end
g(:,1) = [];
g(:,end) = [];

end
