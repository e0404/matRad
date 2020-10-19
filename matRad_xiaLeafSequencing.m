function resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,numOfLevels,visBool)
% multileaf collimator leaf sequencing algorithm 
% for intensity modulated beams with multiple static segments according to 
% Xia et al. (1998) Medical Physics
% 
% call
%   resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,numOfLevels)
%   resultGUI = matRad_xiaLeafSequencing(resultGUI,stf,dij,numOfLevels,visBool)
%
% input
%   resultGUI:          resultGUI struct to which the output data will be added, if
%                       this field is empty resultGUI struct will be created
%   stf:                matRad steering information struct
%   dij:                matRad's dij matrix
%   numOfLevels:        number of stratification levels
%   visBool:            toggle on/off visualization (optional)
%
% output
%   resultGUI:          matRad result struct containing the new dose cube as well as 
%                       the corresponding weights
%
% References
%   [1] http://online.medphys.org/resource/1/mphya6/v25/i8/p1424_s1
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

mode = 'rl'; % sliding window (sw) or reducing level (rl)

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
    wOfCurrBeams = resultGUI.w(1+offset:numOfRaysPerBeam+offset);%REVIEW OFFSET
    
    X = ones(numOfRaysPerBeam,1)*NaN;
    Z = ones(numOfRaysPerBeam,1)*NaN;
        
    for j=1:stf(i).numOfRays
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
    
    % Calculate X and Z position of every fluence's matrix spot
    % z axis = axis of leaf movement!
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
    
    % Save the maximun intensity (Equation 5)
    L_k = max(D_k(:));
    
    % Save the maximun initial intensity matrix value in L_0.
    L_0 = L_k;
    
    % Set k=0, this variable is used for residuals intensity matrices D_k.
    k = 0;
    
    % container to remember generated shapes; allocate space for 10000 shapes
    shapes = NaN*ones(dimOfFluenceMxZ,dimOfFluenceMxX,10000);
    
    if visBool
        clf(seqFig);
        colormap(seqFig,'jet');
        
        seqSubPlots(1) = subplot(2,2,1,'parent',seqFig);
        imagesc(D_k,'parent',seqSubPlots(1));
        set(seqSubPlots(1),'CLim',[0 L_0],'YDir','normal');
        title(seqSubPlots(1),['Beam # ' num2str(i) ': L_0 = ' num2str(L_0) ' - ' num2str(numel(unique(D_0))) ' intensity levels'])
        xlabel(seqSubPlots(1),'x - direction parallel to leaf motion ')
        ylabel(seqSubPlots(1),'z - direction perpendicular to leaf motion ')
        colorbar;
        drawnow
    end
    
    % start sequencer
    while L_k > 0
        
        k = k + 1;
        
        %Plot residual intensity matrix.
        if visBool
            seqSubPlots(2) = subplot(2,2,2,'parent',seqFig);
            imagesc(D_k,'parent',seqSubPlots(2));
            set(seqSubPlots(2),'CLim',[0 L_0],'YDir','normal');
            title(seqSubPlots(2),['k = ' num2str(k) ' - ' num2str(numel(unique(D_k))) ' intensity levels remaining...']);
            xlabel(seqSubPlots(2),'x - direction parallel to leaf motion ');
            ylabel(seqSubPlots(2),'z - direction perpendicular to leaf motion ');
            colorbar
            drawnow
        end
        
        %Rounded off integer. Equation 7.
        m = floor(log2(L_k));
        
        % Convert m=1 if is less than 1. This happens when L_k belong to ]0,2[
        if m < 1
            m = 1;
        end
        
        %Calculate the delivery intensity unit. Equation 6.
        d_k = floor(2^(m-1));
        
        % Opening matrix.
        openingMx = D_k >= d_k;
        
        % Plot opening matrix.
        if visBool
            seqSubPlots(3) = subplot(2,2,3,'parent',seqFig);
            imagesc(openingMx,'parent',seqSubPlots(3));
            set(seqSubPlots(3),'YDir','normal')
            xlabel(seqSubPlots(3),'x - direction parallel to leaf motion ')
            ylabel(seqSubPlots(3),'z - direction perpendicular to leaf motion ')
            title(seqSubPlots(3),'Opening matrix');
            drawnow
        end
        
        if strcmp(mode,'sw') % sliding window technique!
            for j = 1:dimOfFluenceMxZ
                openIx = find(openingMx(j,:) == 1,1,'first');
                if ~isempty(openIx)
                    closeIx = find(openingMx(j,openIx+1:end) == 0,1,'first');
                    if ~isempty(closeIx)
                        openingMx(j,openIx+closeIx:end) = 0;
                    end
                end
                
            end
        elseif strcmp(mode,'rl') % reducing levels technique!
            for j = 1:dimOfFluenceMxZ
                [maxVal,maxIx] = max(openingMx(j,:) .* D_k(j,:));
                if maxVal > 0
                    closeIx = maxIx + find(openingMx(j,maxIx+1:end) == 0,1,'first');
                    if ~isempty(closeIx)
                        openingMx(j,closeIx:end) = 0;
                    end
                    openIx = find(openingMx(j,1:maxIx-1) == 0,1,'last');
                    if ~isempty(openIx)
                        openingMx(j,1:openIx) = 0;
                    end
                end
                                
            end
            
        end
        
        shape_k       = openingMx * d_k;
                
                if visBool
                    seqSubPlots(4) = subplot(2,2,4,'parent',seqFig);
                    imagesc(shape_k,'parent',seqSubPlots(4));
                    set(seqSubPlots(4),'YDir','normal')
                    hold(seqSubPlots(4),'on');
                    xlabel(seqSubPlots(4),'x - direction parallel to leaf motion ')
                    ylabel(seqSubPlots(4),'z - direction perpendicular to leaf motion ')
                    title(seqSubPlots(4),['beam # ' num2str(i) ' shape # ' num2str(k) ' d_k = ' num2str(d_k)]);
                    for j = 1:dimOfFluenceMxZ
                       leftLeafIx = find(shape_k(j,:)>0,1,'first');
                       rightLeafIx = find(shape_k(j,:)>0,1,'last');
                       if leftLeafIx > 1
                           plot(seqSubPlots(4),[.5 leftLeafIx-.5],j-[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),[.5 leftLeafIx-.5],j+[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),[ leftLeafIx-.5 leftLeafIx-.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                       if rightLeafIx<dimOfFluenceMxX
                           plot(seqSubPlots(4),[dimOfFluenceMxX+.5 rightLeafIx+.5],j-[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),[dimOfFluenceMxX+.5 rightLeafIx+.5],j+[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),[ rightLeafIx+.5 rightLeafIx+.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                       if isempty(rightLeafIx) && isempty (leftLeafIx)
                           plot(seqSubPlots(4),[dimOfFluenceMxX+.5 .5],j-[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),[dimOfFluenceMxX+.5 .5],j+[.5 .5] ,'w','LineWidth',2)
                           plot(seqSubPlots(4),0.5*dimOfFluenceMxX*[1 1]+[0.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                    end
                    
                    axis(seqSubPlots(4),'tight');
                    drawnow
                    pause(1);
                end 

        shapes(:,:,k) = shape_k;
        shapesWeight(k) = d_k;
        D_k = D_k - shape_k;
        
        L_k = max(D_k(:)); % eq 5
        
    end
    
    sequencing.beam(i).numOfShapes  = k;
    sequencing.beam(i).shapes       = shapes(:,:,1:k);
    sequencing.beam(i).shapesWeight = shapesWeight(1:k)/numOfLevels*calFac;
    sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
    sequencing.beam(i).fluence      = D_0;
    
    sequencing.w(1+offset:numOfRaysPerBeam+offset,1) = D_0(indInFluenceMx)/numOfLevels*calFac;

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

