function result_Sequencing = matRad_xiaLeafSequencing(w,stf,pln,numOfLevels,visBool)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Multileaf collimator leaf sequencing algorithm for intensity modulated 
% beams with multiple static segments' by Xia and Verhey
% -> http://online.medphys.org/resource/1/mphya6/v25/i8/p1424_s1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if visBool not set toogle off visualization
if nargin < 5
    visBool = 0;
end

mode = 'rl'; % sliding window (sw) or reducing level (rl)

result_Sequencing.w = NaN*w;

numOfBeams = numel(stf);

offset = 0;

for i = 1:numOfBeams
    
    numOfRaysPerBeam = stf(i).numOfRays;
    
    % get relevant weights for current beam
    wOfCurrBeams = w(1+offset:numOfRaysPerBeam+offset);%REVIEW OFFSET
    
    X = zeros(numOfRaysPerBeam,1)*NaN;
    Z = zeros(numOfRaysPerBeam,1)*NaN;
    
    
    for j=1:stf(i).numOfRays
      X(j) = stf(i).ray(j).rayPos_bev(:,1);
      Z(j) = stf(i).ray(j).rayPos_bev(:,3);
    end
    
    
    % sort bixels into matrix
    minX = min(X);
    maxX = max(X);
    minZ = min(Z);
    maxZ = max(Z);
    
    dimOfFluenceMxX = (maxX-minX)/pln.bixelWidth + 1;
    dimOfFluenceMxZ = (maxZ-minZ)/pln.bixelWidth + 1;
    
    %Create the fluence matrix.
    fluenceMx = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    
    % Calculate X and Z positio of every fluence's matrix spot
    % z axis = axis of leaf movement!
    xPos = (X-minX)/pln.bixelWidth+1;
    zPos = (Z-minZ)/pln.bixelWidth+1;
    
    % Make subscripts for fluence matrix
    indInFluenceMx = zPos + (xPos-1)*dimOfFluenceMxZ;
    
    %Save weights in fluence matrix.
    fluenceMx(indInFluenceMx) = wOfCurrBeams;
    
    %imagesc(fluenceMx)
    %xlabel('x - direction parallel to leaf motion ')
    %ylabel('z - direction perpendicular to leaf motion ')
    
    % prepare sequencer
    calFac = max(fluenceMx(:));

    %Stratification
    D_k = round(fluenceMx/calFac*numOfLevels); 
    
    % Save the stratification in the initial intensity matrix D_0.
    D_0 = D_k;
    
    % Save the maximun intensity (Equation 5)
    L_k = max(D_k(:));
    
    % Save the maximun initial intensity matrix value in L_0.
    L_0 = L_k;
    
    % Set k=0, this variable is used for residuals intensity matrices D_k.
    k = 0;
    
    shapes = NaN*ones(dimOfFluenceMxZ,dimOfFluenceMxX,10000);%ANSWER TO MARK WHY 10.000!!!!
    
    if visBool
        clf
        colormap jet
        subplot(2,2,1)
        imagesc(D_k,[0 L_0])
        title(['Beam # ' num2str(i) ': L_0 = ' num2str(L_0) ' - ' num2str(numel(unique(D_0))) ' intensity levels'])
        xlabel('x - direction parallel to leaf motion ')
        ylabel('z - direction perpendicular to leaf motion ')
        colorbar
        drawnow
    end
    
    % start sequencer
    while L_k > 0
        
        k = k + 1;
        
        %Plot residual intensity matrix.
        if visBool
            subplot(2,2,2)
            imagesc(D_k,[0 L_0]);
            title(['k = ' num2str(k) ' - ' num2str(numel(unique(D_k))) ' intensity levels remaining...'])
            xlabel('x - direction parallel to leaf motion ')
            ylabel('z - direction perpendicular to leaf motion ')
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
            subplot(2,2,3)
            imagesc(openingMx);
            xlabel('x - direction parallel to leaf motion ')
            ylabel('z - direction perpendicular to leaf motion ')
            title('Opening matrix');
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
        
                if visBool
                    subplot(2,2,4)
                    imagesc(openingMx)
                    hold on
                    xlabel('x - direction parallel to leaf motion ')
                    ylabel('z - direction perpendicular to leaf motion ')
                    title(['d_k = ' num2str(d_k)]);
                    for j = 1:dimOfFluenceMxZ
                       leftLeafIx = find(openingMx(j,:)>0,1,'first');
                       rightLeafIx = find(openingMx(j,:)>0,1,'last');
                       if leftLeafIx > 1
                           plot([.5 leftLeafIx-.5],j-[.5 .5] ,'w','LineWidth',2)
                           plot([.5 leftLeafIx-.5],j+[.5 .5] ,'w','LineWidth',2)
                           plot([ leftLeafIx-.5 leftLeafIx-.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                       if rightLeafIx<dimOfFluenceMxX
                           plot([dimOfFluenceMxX+.5 rightLeafIx+.5],j-[.5 .5] ,'w','LineWidth',2)
                           plot([dimOfFluenceMxX+.5 rightLeafIx+.5],j+[.5 .5] ,'w','LineWidth',2)
                           plot([ rightLeafIx+.5 rightLeafIx+.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                       if isempty(rightLeafIx) && isempty (leftLeafIx)
                           plot([dimOfFluenceMxX+.5 .5],j-[.5 .5] ,'w','LineWidth',2)
                           plot([dimOfFluenceMxX+.5 .5],j+[.5 .5] ,'w','LineWidth',2)
                           plot(.5*dimOfFluenceMxX*[1 1]+[0.5],j+[.5 -.5] ,'w','LineWidth',2)
                       end
                    end
                    
                    axis tight
                    drawnow
                    pause(1);
                end 

        shape_k       = openingMx * d_k;
        shapes(:,:,k) = shape_k;
        shapesWeight(k) = d_k;
        D_k = D_k - shape_k;
        
        L_k = max(D_k(:)); % eq 5
        
    end
    
    result_Sequencing.beam(i).numOfShapes  = k;
    result_Sequencing.beam(i).shapes       = shapes(:,:,1:k);
    result_Sequencing.beam(i).shapesWeight = shapesWeight(1:k)/numOfLevels*calFac;
    result_Sequencing.beam(i).bixelIx      = 1+offset:numOfRaysPerBeam+offset;
    result_Sequencing.beam(i).fluence      = D_0;
    
    result_Sequencing.w(1+offset:numOfRaysPerBeam+offset) = D_0(indInFluenceMx)/numOfLevels*calFac;

    offset = offset + numOfRaysPerBeam;


end

end

