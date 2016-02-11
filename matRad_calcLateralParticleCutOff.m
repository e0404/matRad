function [ machine ] = matRad_calcLateralParticleCutOff(machine,cutOffLevel,visBool)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate a depth dependend lateral cutoff for each 
% pristine particle beam
% 
% call
%   [ machine ] = matRad_calcLateralParticleCutOff( machine,CutOffLevel,visBool )
%
% input
%   machine:         machine base data file
%   CutOffLevel:     cut off level - number between 0 and 1
%   visBool:         toggle visualization (optional)
%
% output
%   machine:         machine base data file including an additional field representing the lateral
%                    cutoff
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
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

if (cutOffLevel < 0 || cutOffLevel > 0.9999) && (cutOffLevel ~= 1)
   warning('lateral cutoff is out of range - using default cut off of 0.99') 
   cutOffLevel = 0.99;
end

% define some variables needed for the cutoff calculation
vX     = [0 logspace(-1,4,1000)]; % [mm]
%vX     = linspace(0,1000,1000); % [mm]
% integration steps
r_mid   = 0.5*(vX(1:end-1) +  vX(2:end)); % [mm]
dr      = vX(2:end) - vX(1:end-1);
% number of depth points for which a lateral cutoff is determined
NumDepthVal = 30; 

% define function handles for single and double gauss
SG =  @(vR,Sigma)((1/(2*pi*Sigma^2)).*exp(-(vR.^2)./(2*Sigma^2)));
DG =  @(vR,Z,w,Sigma1,Sigma2) Z*(((1-w)*SG(vR,Sigma1)) + (w*SG(vR,Sigma2)));
                         
% loop over all entries in the machine.data struct
for energyIx = 1:length(machine.data)

    % set default depth cut off - finite value will be set during first
    % iteration
    depthDoseCutOff = inf;

    % get indices for which a lateral cutoff should be calculated - always include peak position 
    [~,peakIdx] = max(machine.data(energyIx).Z);
    Idx = round(linspace(1,length(machine.data(energyIx).depths),NumDepthVal-1));
    Idx = unique(sort([Idx peakIdx]));

    for j = 1:length(Idx)
        
        % save depth value
        machine.data(energyIx).LatCutOff.depths(j) = machine.data(energyIx).depths(Idx(j));
        % relative contribution
        relContrib = machine.data(energyIx).Z(Idx(j))/max(machine.data(energyIx).Z);
         
        if strcmp(machine.meta.dataType,'singleGauss')
                    
                    Sigma1 = machine.data(energyIx).sigma(Idx(j));
                    y_r = SG(r_mid,Sigma1);
                 
        elseif strcmp(machine.meta.dataType,'doubleGauss')
            
                    w      = machine.data(energyIx).weight(Idx(j));
                    Sigma1 = machine.data(energyIx).sigma1(Idx(j));
                    Sigma2 = machine.data(energyIx).sigma2(Idx(j));
                    y_r    = DG(r_mid,1,w,Sigma1,Sigma2);
        else
            error('unknown dataType');
        end
        
        if cutOffLevel == 1
            machine.data(energyIx).LatCutOff.CompFac = 1;
            machine.data(energyIx).LatCutOff.CutOff  = Inf;
        else
            % shell integration
            cumArea = cumsum(2*pi.*r_mid.*y_r.*dr);
            
            % check if  -relative contribution in slice is smaller than cut off level
            %           -current depth is behind the peak
            if relContrib < depthDoseCutOff && ...
               machine.data(energyIx).LatCutOff.depths(j) > machine.data(energyIx).peakPos
                  currCutOffLevel = 0;
            else
                  currCutOffLevel = cutOffLevel;
            end
            
            % interpolate cut off
            if cumArea(end) < currCutOffLevel
                error('interpolation impossible. desired cutoff to large for considered area\n');
            end
            
            [cumAreaUnique,IdxUnique] = unique(cumArea);
            try
                if currCutOffLevel == 0
                    r_cut = 0;
                else
                    r_cut = interp1(cumAreaUnique,r_mid(IdxUnique),currCutOffLevel);
                end
            catch
                error('code goes here');
            end
            
%             endIx = find(cumArea == 1,1,'first');
%             if isempty(endIx)
%                 endIx = numel(cumArea);
%             end
%             r_cut = interp1(cumArea(1:endIx),r_mid(1:endIx),currCutOffLevel);
            
            machine.data(energyIx).LatCutOff.CutOff(j) = r_cut;
            
            % set DepthDoseCutOff according to cutoff at entrance dose;
            if j == 1 && strcmp(machine.meta.dataType,'singleGauss')
                depthDoseCutOff = SG(r_cut,Sigma1); 
            elseif j == 1 && strcmp(machine.meta.dataType,'doubleGauss')
                depthDoseCutOff = DG(r_cut,1,w,Sigma1,Sigma2); 
            end
            
            % ensure a monotone increasing lateral cutoff before the bragg peak
            if j > 1 &&  r_cut <  machine.data(energyIx).LatCutOff.CutOff(j-1)
                  machine.data(energyIx).LatCutOff.CutOff(j) =  machine.data(energyIx).LatCutOff.CutOff(j-1);
            end

            % Compenstaion factor to rescale the dose within the cut off in order not to lose integral dose
            machine.data(energyIx).LatCutOff.CompFac = 1/cutOffLevel;
            
        end
   end
    
end    



            
%% visualization

if visBool
    
    % get random energy index
    energyIx = round((length(machine.data)-1).*rand(1,1) + 1);
    % set depth position - 1 means plotting the entry profile
    j = Idx(1);
    CutOff = machine.data(energyIx).LatCutOff.CutOff(j);
    if isinf(CutOff)
        return
    end
    % plot 3D cutoff at one specific depth on a rough grid
    Step = 0.5;
    vLatX = -CutOff*1.5 : Step : CutOff*1.5; % [mm]
    midPos = round(length(vLatX)/2);
    [X,Y] = meshgrid(vLatX,vLatX);
    vRadDist = sqrt(X.^2 + Y.^2);
    
   
    if strcmp(machine.meta.dataType,'singleGauss')
            vDose = SG(vRadDist,machine.data(energyIx).sigma(j));    
    elseif strcmp(machine.meta.dataType,'doubleGauss')
            vDose = DG(vRadDist,1,...
             machine.data(energyIx).weight(j),...
             machine.data(energyIx).sigma1(j),...
             machine.data(energyIx).sigma2(j));
    end
         
    [~,LevelIdx] = min(abs(X(1,:)-CutOff));
    
    DoseLevel = vDose(midPos,LevelIdx);
    
    figure,set(gcf,'Color',[1 1 1]);
    subplot(221),surf(X,Y,vDose),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'), hold on, axis tight
    contour3(X,Y,vDose,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r'),hold on;
    title({['beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)], ['cutoff = ' num2str(cutOffLevel)]}),set(gca,'FontSize',12);
    
    subplot(222),surf(X,Y,vDose),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'), hold on, axis tight
    contour3(X,Y,vDose,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r'),hold on; title(['intensity profile; cutoff = ' num2str(cutOffLevel)]),view(0,90)
    
    if strcmp(machine.meta.dataType,'singleGauss')
         vDoseLat =  machine.data(energyIx).Z(j)*SG(vLatX,machine.data(energyIx).sigma(j));
    elseif strcmp(machine.meta.dataType,'doubleGauss')
         vDoseLat =  DG(vLatX,machine.data(energyIx).Z(j),machine.data(energyIx).weight(j),...
         machine.data(energyIx).sigma1(j),machine.data(energyIx).sigma2(j));
    end
    
    subplot(223),plot(vLatX,vDoseLat,'LineWidth',3),grid on, grid minor, hold on
    plot([CutOff,CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on
    plot([-CutOff,-CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on, title('lateral profile 2D - cross section')
    
    subplot(224),surf(X,Y,vDose),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'),colormap(parula(256)), hold on
    title(['proton beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)]),set(gca,'FontSize',12);
    contour3(X,Y,vDose,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r');title('lateral profile 3D'); view([0 0]);

    % generate 4 depth points for visualization
    idx = round(linspace(1,length(machine.data(energyIx).LatCutOff.depths),4));
    vDepth = machine.data(energyIx).depths;
  
    figure,set(gcf,'Color',[1 1 1]);
    
    maxZ = max(machine.data(energyIx).Z);
    subplot(321),plot(vDepth, machine.data(energyIx).Z,'LineWidth',3),
    grid on, grid minor, title(['ddd with energy ' num2str(machine.data(energyIx).energy)])
    hold on,
    plot([machine.data(energyIx).LatCutOff.depths(idx(1)),machine.data(energyIx).LatCutOff.depths(idx(1))],...
        [0 maxZ],'r','LineWidth',2),hold on    
    plot([machine.data(energyIx).LatCutOff.depths(idx(2)),machine.data(energyIx).LatCutOff.depths(idx(2))],...
        [0 maxZ],'k','LineWidth',2),
    plot([machine.data(energyIx).LatCutOff.depths(idx(3)),machine.data(energyIx).LatCutOff.depths(idx(3))],...
        [0 maxZ],'b','LineWidth',2),
    plot([machine.data(energyIx).LatCutOff.depths(idx(4)),machine.data(energyIx).LatCutOff.depths(idx(4))],...
        [0 maxZ],'g','LineWidth',2),
    set(gca,'FontSize',12)
    
    subplot(322),
    maxLateralCutOff = max(machine.data(energyIx).LatCutOff.CutOff);
    plot(machine.data(energyIx).LatCutOff.depths, machine.data(energyIx).LatCutOff.CutOff,'LineWidth',3),grid on, grid minor,
    ylabel('lateral cutoff in [mm]')
    title(['desired CutOff is ' num2str(cutOffLevel*100) '%']),hold on
    plot([machine.data(energyIx).LatCutOff.depths(idx(1)),machine.data(energyIx).LatCutOff.depths(idx(1))],...
        [0 maxLateralCutOff],'r','LineWidth',2),hold on    
    plot([machine.data(energyIx).LatCutOff.depths(idx(2)),machine.data(energyIx).LatCutOff.depths(idx(2))],...
        [0 maxLateralCutOff],'k','LineWidth',2),
    plot([machine.data(energyIx).LatCutOff.depths(idx(3)),machine.data(energyIx).LatCutOff.depths(idx(3))],...
        [0 maxLateralCutOff],'b','LineWidth',2),
    plot([machine.data(energyIx).LatCutOff.depths(idx(4)),machine.data(energyIx).LatCutOff.depths(idx(4))],...
        [0 maxLateralCutOff],'g','LineWidth',2),
    set(gca,'FontSize',12)
    
    % generate fine grid
    Step = 0.05;
    vLatX = -150:Step:150; % [mm]
    midPos = round(length(vLatX)/2);
    [X,Y] = meshgrid(vLatX,vLatX);
    vRadDist = sqrt(X.^2 + Y.^2);
    
    subplot(323),
    [~,refIdx] = min(abs(machine.data(energyIx).depths - machine.data(energyIx).LatCutOff.depths(idx(1))));
    title('intensity profile')
    if strcmp(machine.meta.dataType,'singleGauss')
       vDose = SG(vRadDist,machine.data(energyIx).sigma(refIdx));
    elseif strcmp(machine.meta.dataType,'doubleGauss')
       vDose = DG(vRadDist,1,machine.data(energyIx).weight(refIdx),...
       machine.data(energyIx).sigma1(refIdx),machine.data(energyIx).sigma2(refIdx));
    end
    plot(vLatX,vDose(midPos,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.CutOff(idx(1));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadDist < cutOff )) * Step^2;
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100) '% at depth: ' num2str(machine.data(energyIx).depths(refIdx)) ' mm']),set(gca,'FontSize',12)
    
    subplot(324),
    [~,refIdx] = min(abs(machine.data(energyIx).depths - machine.data(energyIx).LatCutOff.depths(idx(2))));
    if strcmp(machine.meta.dataType,'singleGauss')
          vDose = SG(vRadDist,machine.data(energyIx).sigma(refIdx));
    elseif strcmp(machine.meta.dataType,'doubleGauss')
        vDose = DG(vRadDist,1,machine.data(energyIx).weight(refIdx),...
         machine.data(energyIx).sigma1(refIdx),machine.data(energyIx).sigma2(refIdx));
    end
    plot(vLatX,vDose(midPos,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.CutOff(idx(2));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadDist<cutOff))* Step^2;
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100) '% at depth :' num2str(machine.data(energyIx).depths(refIdx)) ' mm']),set(gca,'FontSize',12)
   
    subplot(325),
    [~,refIdx] = min(abs(machine.data(energyIx).depths - machine.data(energyIx).LatCutOff.depths(idx(3))));
    if strcmp(machine.meta.dataType,'singleGauss')
         vDose = SG(vRadDist,machine.data(energyIx).sigma(refIdx));
     elseif strcmp(machine.meta.dataType,'doubleGauss')
         vDose =  DG(vRadDist,1,machine.data(energyIx).weight(refIdx),...
           machine.data(energyIx).sigma1(refIdx),machine.data(energyIx).sigma2(refIdx));
    end
    plot(vLatX,vDose(midPos,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.CutOff(idx(3));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadDist<cutOff))* Step^2;
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100) '% at depth: ' num2str(machine.data(energyIx).depths(refIdx)) 'mm']),set(gca,'FontSize',12)
 
    
    subplot(326),
    [~,refIdx] = min(abs(machine.data(energyIx).depths - machine.data(energyIx).LatCutOff.depths(idx(4))));
    if strcmp(machine.meta.dataType,'singleGauss')
         vDose = SG(vRadDist,machine.data(energyIx).sigma(refIdx));
    elseif strcmp(machine.meta.dataType,'doubleGauss')
         vDose =  DG(vRadDist,1,machine.data(energyIx).weight(refIdx),...
             machine.data(energyIx).sigma1(refIdx),machine.data(energyIx).sigma2(refIdx));
    end
    plot(vLatX,vDose(midPos,:),'LineWidth',3),grid on, grid minor, hold on
    cutOff = machine.data(energyIx).LatCutOff.CutOff(idx(4));
    plot([-cutOff -cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    plot([cutOff cutOff],[0 max(vDose(:))],'r','LineWidth',2)
    Atot = sum(vDose(vRadDist<cutOff))* Step^2;
    title(['lateral profile;  calc. cutOff: ' num2str((Atot)*100) '% at depth: ' num2str(machine.data(energyIx).depths(refIdx)) ' mm']),set(gca,'FontSize',12)
    
    % plot cutoff of different energies
    figure,set(gcf,'Color',[1 1 1]);
    NumPlots = 20;
    ix = round(linspace(1,length(machine.data),NumPlots));
    for i = 1:20
        plot(machine.data(ix(i)).LatCutOff.depths,machine.data(ix(i)).LatCutOff.CutOff,'LineWidth',1.5),hold on
        cellLegend{i} = [num2str(machine.data(ix(i)).energy) ' MeV'];
    end
    grid on, grid minor,xlabel('depth in [mm]'),ylabel('lateral cutoff in [mm]')
    title(['cutoff level = ' num2str(cutOffLevel)])
    legend(cellLegend)
end



end

