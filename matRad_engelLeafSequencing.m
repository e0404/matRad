function resultGUI = matRad_engelLeafSequencing(resultGUI,stf,dij,numOfLevels,visBool)
% multileaf collimator leaf sequencing algorithm for intensity modulated 
% beams with multiple static segments accroding to Engel et al. 2005
% Discrete Applied Mathematics
% 
% call
%   resultSequencing = matRad_engelSequencing(w,stf,pln,numOfLevels,visBool)
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
%   resultGUI:          matRad result struct containing the new dose cube
%                       as well as the corresponding weights
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0166218X05001411
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
    wOfCurrBeams = resultGUI.w(1+offset:numOfRaysPerBeam+offset);
    
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
    
    % container to remember generated shapes; allocate space for 10000 shapes
    shapes = NaN*ones(dimOfFluenceMxZ,dimOfFluenceMxX,10000);
  
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
    
    % start sequencer
    while max(D_k(:) > 0)
        
        %calculate the difference matrix diffMat
        diffMat = diff([zeros(size(D_k,1),1) D_k zeros(size(D_k,1),1)],[],2);
        
        %calculate complexities
        c = sum(max(0,diffMat),2); %TNMU-row-complexity
        com = max(c); %TNMU complexity
        g = com - c; %row complexity gap
        
        %initialize segment
        segment = zeros(size(D_k));
        
        k = k + 1;
        
        %Plot residual intensity matrix.
        if visBool
            seqSubPlots(2) = subplot(2,2,2,'parent',seqFig);
            imagesc(D_k,'parent',seqSubPlots(2));
            set(seqSubPlots(2),'CLim',[0 numOfLevels],'YDir','normal');
            title(seqSubPlots(2),['k = ' num2str(k)]);
            colorbar
            drawnow
        end
       
        
        %loop over all rows 
        for j=1:size(D_0,1)
            
            %determine essential intervals
            data(j).left(1) = 0; %left interval limit, actual for an empty interval
            data(j).right(1) = 0; %right interal limit, actual for an empty interval
            data(j).v(1) = g(j);  %greatest number such that the inequalities (6) resp. (7) is satisfied with u=v
            data(j).w(1) = inf; %smallest number in the interval
            data(j).u(1) = data(j).v(1); %min(v,w)
            
            [~, pos, ~] = find(diffMat(j,:) > 0); % indices of all positive elements in the j. row of diffmat
            [~, neg, ~] = find(diffMat(j,:) < 0); % indices of all negative elements in the j. row of diffMat
            
            n=2;
            
            %loop over the positive elements in the j. row of diffmat ->
            %possible left interval limits
            for m=1:size(pos,2)
                
                %loop over the negative elements in the j. row of diffMat ->
                %possible right interval limit
                for l=1:size(neg,2)
                    
                    %take only intervals I=[l,r] with l<=r
                    if pos(m) <= neg(l)-1
                        
                        %set interval limits
                        data(j).left(n) = pos(m);
                        data(j).right(n) = neg(l)-1;
                        
                        %calculate v according to Lemma 8
                        if g(j) <= abs( diffMat(j,pos(m)) + diffMat(j,neg(l)) )
                            data(j).v(n) = min( diffMat(j,pos(m)), -diffMat(j,neg(l)) ) + g(j);
                        else
                            data(j).v(n) = ( diffMat(j, pos(m)) - diffMat(j, neg(l)) + g(j)) / 2;
                        end
                        
                        %calculate w and u according to equality (11) and
                        %(12) 
                        data(j).w(n) = min(D_k(j,pos(m):(neg(l)-1)));
                        data(j).u(n) = min(data(j).v(n), data(j).w(n));
                        
                        n = n+1;
                    end
                end
            end
            
            u(j) = max(data(j).u);
            
        end
        
        %calculate u_max from theorem 9
        d_k = min(u);
        
        %loop over all rows
        for j=1:size(D_0,1)
 
            %find all possible (and essential) intervals
            candidate = find(data(j).u >= d_k); 
            
            %calculate the potential of the possible intervals
            
            %initialize p as -Inf
            data(j).p(1:length(data(j).left)) = -Inf;
            
            %loop over all possible intervals
            for s=1:size(candidate,2)
                
                if (s==1 && data(j).left(candidate(s)) == 0)
                    data(j).p(candidate(1)) = 0;
                    
                    
                else
                    %calculate p1 according to equality (17)
                    if (d_k == diffMat(j, data(j).left(candidate(s))) && d_k ~= D_k(j, data(j).left(candidate(s))))
                        p1 = 1;
               
                    else
                        p1 = 0;
                    
                    end
                    
                    %calculate p2 according to equalitiy (18)
                   % if data(j).right(candidate(s)) < size(D_0, 2)
                        
                        if (d_k == -diffMat(j, data(j).right(candidate(s))+1) && d_k ~= D_k(j, data(j).right(candidate(s))))
                            p2 = 1;
                        else
                            p2 = 0;
                        end
                        
%                     else
%                         
%                         if d_k == -diffMat(j, data(j).right(candidate(s))+1)
%                             p2 = 1;
%                         else
%                             p2 = 0;
%                         end
%                         
%                     end                                     
                    
                    %calculate p3 according to equality (19)
                    p3 = size(find(D_k(j, data(j).left(candidate(s)):data(j).right(candidate(s))) == d_k),2);
                    
                    data(j).p(candidate(s)) = p1 + p2+ p3;
                    
                end
                
            end
                
                %determinate intervals with maximum potential
                maxPot = find(data(j).p == max(data(j).p));
                
                %if several intervals have maximum potential, select
                %the interval which has maximum length
                if size(maxPot,2) > 1

                    for t=1:size(maxPot,2)
                        if t==1 && data(j).left(maxPot(t)) == 0
                            data(j).l(1) = 0;
                        else
                            data(j).l(maxPot(t)) = data(j).right(maxPot(t)) - data(j).left(maxPot(t)) + 1;
                        end
                    end
                    
                    %data(j).l(maxPot) = data(j).right(maxPot) - data(j).left(maxPot) + 1;
                    
                    maxLength = find(data(j).l == max(data(j).l));
                     
                    %left and right interval limits of the selected
                    %interval
                    leftIntLimit(j) = data(j).left(maxLength(1));
                    rightIntLimit(j) = data(j).right(maxLength(1));
                    
 
                else
                    
                    %left and right interval limits of the selected
                    %interval
                    leftIntLimit(j) = data(j).left(maxPot);
                    rightIntLimit(j) = data(j).right(maxPot);
                    
                    
                end
                
                %create segment associated by the selected interval
                if leftIntLimit(j) ~= 0 
                    
                    segment(j,leftIntLimit(j):rightIntLimit(j)) = 1;
 
                end

        end
        
        %write the segment in shape_k
        shape_k = segment;

        %show the leaf positions
        if visBool
            seqSubPlots(4) = subplot(2,2,3.5,'parent',seqFig);
            imagesc(shape_k,'parent',seqSubPlots(4));
            hold(seqSubPlots(4),'on');
            set(seqSubPlots(4),'YDir','normal')
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
                    plot(seqSubPlots(4),.5*dimOfFluenceMxX*[1 1]+[0.5],j+[.5 -.5] ,'w','LineWidth',2)
                end
            end
            
            axis tight
            drawnow
            pause(1);
        end
    

        %save shape_k in container
        shapes(:,:,k) =  shape_k;
        
        %save the calculated MU
        shapesWeight(k) = d_k; 
        
        %calculate  new matrix, the  diference matrix and complexities
        D_k = D_k - d_k*shape_k; 
        
        %delete variables
        clear data;
        clear segment;
        clear u;
        clear leftIntLimit;
        clear rightIntLimit;
       
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

