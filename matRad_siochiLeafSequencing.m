function resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,numOfLevels,visBool)
% multileaf collimator leaf sequencing algorithm for intensity modulated
% beams with multiple static segments according to Siochi (1999)
% International Journal of Radiation Oncology * Biology * Physics,
% originally implemented in PLUNC (https://sites.google.com/site/planunc/)
%
% Implemented in matRad by Eric Christiansen, Emily Heath, and Tong Xu
%
% call
%   resultSequencing =
%   matRad_siochiLeafSequencing(w,stf,pln,numOfLevels,visBool)
%
% input
%   resultGUI:          resultGUI struct to which the output data will be
%                       added, if this field is empty resultGUI struct will
%                       be created
%   stf:                matRad steering information struct
%   dij:                matRad's dij matrix
%   numOfLevels:        number of stratification levels
%   visBool:            toggle on/off visualization (optional)
%
% output
%   resultGUI:          matRad result struct containing the new dose cube
%                       as well as the corresponding weights
%
% References
%   [1] https://www.ncbi.nlm.nih.gov/pubmed/10078655
%
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
    
    % get relevant weights for current beam
    wOfCurrBeams = resultGUI.wUnsequenced(1+offset:numOfRaysPerBeam+offset);%REVIEW OFFSET
    
    X = ones(numOfRaysPerBeam,1)*NaN;
    Z = ones(numOfRaysPerBeam,1)*NaN;
    
    for j = 1:stf(i).numOfRays
        X(j) = stf(i).ray(j).rayPos_bev(:,1);
        Z(j) = stf(i).ray(j).rayPos_bev(:,3);
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
    
    % Stratification
    calFac = max(fluenceMx(:));
    D_k = round(fluenceMx/calFac*numOfLevels);
    
    % Save the stratification in the initial intensity matrix D_0.
    D_0 = D_k;
    
    % container to remember generated shapes; allocate space for 10000
    % shapes
    shapes = NaN*ones(dimOfFluenceMxZ,dimOfFluenceMxX,10000);
    shapesWeight = zeros(10000,1);
    k = 0;
    
    if visBool
        clf(seqFig);
        colormap(seqFig,'jet');
        
        seqSubPlots(1) = subplot(2,2,1,'parent',seqFig);
        imagesc(D_k,'parent',seqSubPlots(1));
        set(seqSubPlots(1),'CLim',[0 numOfLevels],'YDir','normal');
        title(seqSubPlots(1),['Beam # ' num2str(i) ': max(D_0) = ' num2str(max(D_0(:))) ' - ' num2str(numel(unique(D_0))) ' intensity levels']);
        xlabel(seqSubPlots(1),'x - direction parallel to leaf motion ')
        ylabel(seqSubPlots(1),'z - direction perpendicular to leaf motion ')
        colorbar;
        drawnow
    end
    
    D_k_nonZero = (D_k~=0);
    [D_k_Z, D_k_X] = ind2sub([dimOfFluenceMxZ,dimOfFluenceMxX],find(D_k_nonZero));
    D_k_MinZ = min(D_k_Z);
    D_k_MaxZ = max(D_k_Z);
    D_k_MinX = min(D_k_X);
    D_k_MaxX = max(D_k_X);
    
    %Decompose the port, do rod pushing
    [tops, bases] = matRad_siochiDecomposePort(D_k,dimOfFluenceMxZ,dimOfFluenceMxX,D_k_MinZ,D_k_MaxZ,D_k_MinX,D_k_MaxX);
    %Form segments with and without visualization
    if visBool
        [shapes,shapesWeight,k,D_k]=matRad_siochiConvertToSegments(shapes,shapesWeight,k,tops,bases,visBool,i,D_k,numOfLevels,seqFig,seqSubPlots);
    else
        [shapes,shapesWeight,k]=matRad_siochiConvertToSegments(shapes,shapesWeight,k,tops,bases);
    end
    
    sequencing.beam(i).numOfShapes  = k;
    sequencing.beam(i).shapes       = shapes(:,:,1:k);
    sequencing.beam(i).shapesWeight = shapesWeight(1:k)/numOfLevels*calFac;
    sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
    sequencing.beam(i).fluence      = D_0;
    sequencing.beam(i).sum          = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    
    for j = 1:k
        sequencing.beam(i).sum = sequencing.beam(i).sum+sequencing.beam(i).shapes(:,:,j)*sequencing.beam(i).shapesWeight(j);
    end
    sequencing.w(1+offset:numOfRaysPerBeam+offset,1) = sequencing.beam(i).sum(indInFluenceMx);
    
    offset = offset + numOfRaysPerBeam;

end

resultGUI.w          = sequencing.w;
resultGUI.wSequenced = sequencing.w;

resultGUI.sequencing   = sequencing;
resultGUI.apertureInfo = matRad_sequencing2ApertureInfo(sequencing,stf);

doseSequencedDoseGrid = reshape(dij.physicalDose{1} * sequencing.w,dij.doseGrid.dimensions);
% interpolate to ct grid for visualiation & analysis
resultGUI.physicalDose = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                        doseSequencedDoseGrid, ...
                                        dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z);

% if weights exists from an former DAO remove it
if isfield(resultGUI,'wDao')
    resultGUI = rmfield(resultGUI,'wDao');
end

end

function [tops, bases] = matRad_siochiDecomposePort(map,dimZ,dimX,minZ,maxZ,minX,maxX)
%Returns tops and bases of a fluence matrix "map" for Siochi leaf
%sequencing algorithm (rod pushing part).  Accounts for collisions and
%tongue and groove (Tng) effects.

tops = zeros(dimZ, dimX);
bases = zeros(dimZ, dimX);

for i = minX:maxX
    maxTop = -1;
    TnG = 1;
    for j = minZ:maxZ
        if i == minX
            bases(j,i) = 1;
            tops(j,i) = bases(j,i)+map(j,i)-1;
        else %assign trial base positions
            if map(j,i) >= map(j,i-1) %current rod >= previous, match the bases
                bases(j,i) = bases(j,i-1);
                tops(j,i) = bases(j,i)+map(j,i)-1;
            else %current rod <previous
                if map(j,i) == 0 %rod length=0, put in in next slab after top of previous
                    bases(j,i) = tops(j,i-1)+1;
                    tops(j,i) = bases(j,i)-1;
                else %rod length~=0, match tops
                    tops(j,i) = tops(j,i-1);
                    bases(j,i) = tops(j,i)-map(j,i)+1;
                end
            end
        end
        %determine which rod has the highest top in column
        if tops(j,i) > maxTop
            maxTop = tops(j,i);
            maxRow = j;
        end
    end
    
    %Correct for collision and tongue and groove error
    while(TnG)
        %go from maxRow down checking for TnG.  This occurs when a shorter
        %rod is "peeking over" a longer one in the direction transverse to
        %the leaf motion.  To fix this, match either the tops or bases of
        %the rods.
        for j = (maxRow-1):-1:minZ
            if map(j,i) < map(j+1,i)
                if tops(j,i) > tops(j+1,i)
                    tops(j+1,i) = tops(j,i);
                    bases(j+1,i) = tops(j+1,i)-map(j+1,i)+1;
                elseif bases(j,i) < bases(j+1,i)
                    bases(j,i) = bases(j+1,i);
                    tops(j,i) = bases(j,i)+map(j,i)-1;
                end
            else
                if tops(j,i) < tops(j+1,i)
                    tops(j,i) = tops(j+1,i);
                    bases(j,i) = tops(j,i)-map(j,i)+1;
                elseif bases(j,i) > bases(j+1,i)
                    bases(j+1,i) = bases(j,i);
                    tops(j+1,i) = bases(j+1,i)+map(j+1,i)-1;
                end
            end
        end
        %go from maxRow up checking for TnG
        for j = (maxRow+1):maxZ
            if map(j,i) < map(j-1,i)
                if tops(j,i) > tops(j-1,i)
                    tops(j-1,i) = tops(j,i);
                    bases(j-1,i) = tops(j-1,i)-map(j-1,i)+1;
                elseif bases(j,i) < bases(j-1,i)
                    bases(j,i) = bases(j-1,i);
                    tops(j,i) = bases(j,i)+map(j,i)-1;
                end
            else
                if tops(j,i) < tops(j-1,i)
                    tops(j,i) = tops(j-1,i);
                    bases(j,i) = tops(j,i)-map(j,i)+1;
                elseif bases(j,i) > bases(j-1,i)
                    bases(j-1,i) = bases(j,i);
                    tops(j-1,i) = bases(j-1,i)+map(j-1,i)-1;
                end
            end
        end
        %now check if all TnG conditions have been removed
        TnG = 0;
        for j = (minZ+1):maxZ
            if map(j,i) < map(j-1,i);
                if tops(j,i) > tops(j-1,i)
                    TnG = 1;
                elseif bases(j,i) < bases(j-1,i)
                    TnG = 1;
                end
            else
                if tops(j,i) < tops(j-1,i)
                    TnG = 1;
                elseif bases(j,i) > bases(j-1,i)
                    TnG = 1;
                end
            end
        end
    end
end

end

function [shapes,shapesWeight,k,D_k] = matRad_siochiConvertToSegments(shapes,shapesWeight,k,tops,bases,visBool,i,D_k,numOfLevels,seqFig,seqSubPlots)
%Convert tops and bases to shape matrices.  These are taken as to be the
%shapes of uniform level/elevation after the rods are pushed.
if nargin < 6
    visBool = 0;
end


levels = max(tops(:));

for level = 1:levels
    %check if slab is new
    if matRad_siochiDifferentSlab(tops,bases,level)
        k = k+1; %increment number of unique slabs
        shape_k = (bases <= level).*(level <= tops); %shape of current slab
        shapes(:,:,k) = shape_k;
    end
    shapesWeight(k) = shapesWeight(k)+1; %if slab is not unique, this increments weight again
    
    if visBool
        %show the leaf positions
        [dimZ,dimX] = size(tops);
        seqSubPlots(4) = subplot(2,2,3.5,'parent',seqFig);
        imagesc(shape_k,'parent',seqSubPlots(4));
        hold(seqSubPlots(4),'on');
        set(seqSubPlots(4),'YDir','normal')
        xlabel(seqSubPlots(4),'x - direction parallel to leaf motion ')
        ylabel(seqSubPlots(4),'z - direction perpendicular to leaf motion ')
        title(seqSubPlots(4),['beam # ' num2str(i) ' shape # ' num2str(k) ' d_k = ' num2str(shapesWeight(k))]);
        for j = 1:dimZ
            leftLeafIx = find(shape_k(j,:)>0,1,'first');
            rightLeafIx = find(shape_k(j,:)>0,1,'last');
            if leftLeafIx > 1
                plot(seqSubPlots(4),[.5 leftLeafIx-.5],j-[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),[.5 leftLeafIx-.5],j+[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),[ leftLeafIx-.5 leftLeafIx-.5],j+[.5 -.5] ,'w','LineWidth',2)
            end
            if rightLeafIx<dimX
                plot(seqSubPlots(4),[dimX+.5 rightLeafIx+.5],j-[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),[dimX+.5 rightLeafIx+.5],j+[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),[ rightLeafIx+.5 rightLeafIx+.5],j+[.5 -.5] ,'w','LineWidth',2)
            end
            if isempty(rightLeafIx) && isempty (leftLeafIx)
                plot(seqSubPlots(4),[dimX+.5 .5],j-[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),[dimX+.5 .5],j+[.5 .5] ,'w','LineWidth',2)
                plot(seqSubPlots(4),.5*dimX*[1 1]+[0.5],j+[.5 -.5] ,'w','LineWidth',2)
            end
        end
        pause(1);
        
        %Plot residual intensity matrix.
        D_k = D_k-shape_k; %residual intensity matrix for visualization
        seqSubPlots(2) = subplot(2,2,2,'parent',seqFig);
        imagesc(D_k,'parent',seqSubPlots(2));
        set(seqSubPlots(2),'CLim',[0 numOfLevels],'YDir','normal');
        title(seqSubPlots(2),['k = ' num2str(k)]);
        colorbar
        drawnow
        
        axis tight
        drawnow
    end
end

end

function diffSlab = matRad_siochiDifferentSlab(tops,bases,level)
%Returns 1 if slab level is different than slab level-1 0 otherwise

if level == 1 %first slab is automatically different
    diffSlab = 1;
else
    shapeLevel = (bases <= level).*(level <= tops); %shape of slab with current level
    shapeLevel_1 = (bases <= level-1).*(level-1 <= tops); %shape of slab with previous level
    diffSlab = ~isequal(shapeLevel,shapeLevel_1); %tests if slabs are equal; isequaln was not giving correct results
end

end

