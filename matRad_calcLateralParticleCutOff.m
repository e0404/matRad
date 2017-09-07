function [ machine ] = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf,visBool)
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
TypeOfCutOffCalc = '2D'; % 'classic','1D','2D','3D',
LcutDefaultSigma = 3.5;  % in sigma units - this value is only used by the classica approach - only 0.2% of fluence is missed when using 3.5
conversionFactor = 1.6021766208e-02;

% function handle for calculating depth dose
SumGauss = @(x,mu,SqSigma,w) ((1./sqrt(2*pi*ones(numel(x),1) * SqSigma') .* ...
                              exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) * SqSigma' ))) * w);
                          
if (cutOffLevel < 0 || cutOffLevel > 0.9999) && (cutOffLevel ~= 1)
   warning('lateral cutoff is out of range - using default cut off of 0.99') 
   cutOffLevel = 0.99;
end
% define some variables needed for the cutoff calculation
vX     = [0 logspace(-1,5,1300)]; % [mm]

% integration steps
r_mid   = (0.5*(vX(1:end-1) +  vX(2:end)))'; % [mm]
dr      = (vX(2:end) - vX(1:end-1))';

% number of depth points for which a lateral cutoff is determined
NumDepthVal = 30; 

% define function handles for single and double gauss
SG    =  @(vR,Sigma)((1/(2*pi*Sigma^2)).*exp(-(vR.^2)./(2*Sigma^2)));
DG    =  @(vR,Z,w,Sigma1,Sigma2) Z*(((1-w)*SG(vR,Sigma1)) + (w*SG(vR,Sigma2)));
CF    =  @(LcutSigma)(1/(1-exp((-LcutSigma^2)/2)));
CDF   =  @(vR,Si) (1/2)*(1+erf(vR/(Si*(sqrt(2)))));

% extract SSD for each bixel
vSSDBixel = ones(1,length([stf.ray(:).energy]));
cnt = 1;
for i  = 1:length(stf.ray)
    vSSDBixel(cnt:cnt+numel([stf.ray(i).energy])-1) = stf.ray(i).SSD;
    cnt = cnt + numel(stf.ray(i).energy);
end

% setup energy sigma look up table
energySigmaLUT  = unique([[stf.ray(:).energy]; [stf.ray(:).focusIx] ; vSSDBixel]','rows');
rangeShifterLUT = [stf.ray(:).rangeShifter];

% calculate for each energy its inital beam width considering foci and SSD
for i = 1:size(energySigmaLUT,1)
    energyIx = find(ismember([machine.data(:).energy],energySigmaLUT(i,1)));
    energySigmaLUT(i,4) = matRad_interp1(machine.data(energyIx).initFocus.dist(energySigmaLUT(i,2),:)',...
                                  machine.data(energyIx).initFocus.sigma(energySigmaLUT(i,2),:)',...
                                  energySigmaLUT(i,3));
end

% find for each energy the broadest inital beam width
uniqueEnergies = unique(energySigmaLUT(:,1));
largestFocus4uniqueEnergies = NaN * ones(numel(uniqueEnergies),1);
ix_Max                      = NaN * ones(numel(uniqueEnergies),1);
for i = 1:numel(uniqueEnergies)
    [largestFocus4uniqueEnergies(i), ix_Max(i)] = max(energySigmaLUT(uniqueEnergies(i) == energySigmaLUT(:,1),4));
end

% get energy indices for looping
vEnergiesIx = find(ismember([machine.data(:).energy],uniqueEnergies(:,1)));

cnt = 0;    

% loop over all entries in the machine.data struct
for energyIx = vEnergiesIx
   
    % set default depth cut off - finite value will be set during first
    % iteration
    depthDoseCutOff = inf;

    % get indices for which a lateral cutoff should be calculated - always include peak position 
    if isstruct(machine.data(energyIx).Z)
        idd = SumGauss(machine.data(energyIx).depths,machine.data(energyIx).Z.mean,...
                                                     machine.data(energyIx).Z.width.^2,...
                                                     machine.data(energyIx).Z.weight);
    else
        idd = machine.data(energyIx).Z;
    end
    
    [~,peakixDepth] = max(idd); 
    % define depth positions for which lateral cut off should be calculated
    ixDepth = round(linspace(1,length(machine.data(energyIx).depths),NumDepthVal-1));
    
    
    ixDepth = unique(sort([ixDepth peakixDepth]));
    
    % get inital beam width
    cnt = cnt +1 ;
    % % calculate maximum dose in spot
    baseData                   = machine.data(energyIx);
    radDepths                  = machine.data(energyIx).depths(peakixDepth) ;
    radialDist_sq              = 0;
    maxfocusIx                 = energySigmaLUT(ix_Max(cnt),2);
    maxSSD                     = energySigmaLUT(ix_Max(cnt),3);
    radiationMode              = stf(1).radiationMode;
    rangeShifter               = rangeShifterLUT(ix_Max(cnt));
    baseData.LatCutOff.CompFac = 1;
    
    dosePeakPos = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, maxSSD, maxfocusIx, baseData, rangeShifter, radiationMode);   
    iddPeak     = idd(peakixDepth) * conversionFactor;    

     radialDist_sq  = r_mid.^2;
     
    for j = 1:length(ixDepth)
        
        % save depth value
        machine.data(energyIx).LatCutOff.depths(j) = machine.data(energyIx).depths(ixDepth(j));
        radDepths      = machine.data(energyIx).LatCutOff.depths(j) * ones(numel(r_mid),1);       
        dose_r         = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, maxSSD, maxfocusIx, baseData, rangeShifter, radiationMode);
       
        if cutOffLevel == 1
            machine.data(energyIx).LatCutOff.CompFac = 1;
            machine.data(energyIx).LatCutOff.CutOff  = Inf;
        else
            
            switch TypeOfCutOffCalc

                case 'classic'
                      if strcmp(machine.meta.dataType,'singleGauss')
                        LcutDefaultSigma1 = LcutDefaultSigma;
                        IX = find(r_mid  >  Sigma1 * LcutDefaultSigma1 ,1 ,'first'); 
                        machine.data(energyIx).LatCutOff.CompFac = CF(LcutDefaultSigma1);
                      elseif strcmp(machine.meta.dataType,'doubleGauss')
                        LcutDefaultSigma2 = 0.85*LcutDefaultSigma2;
                        IX   =  find(r_mid  >  Sigma2 * 0.85*LcutDefaultSigma2 ,1 ,'first');   
                        machine.data(energyIx).LatCutOff.CompFac = CF(LcutDefaultSigma2);
                      end  

                case '1D'
                    
                      IX = find(dose_r  <= (1-cutOffLevel) * dosePeakPos ,1 ,'first');   
                      machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;
                      
                case '2D'
                      cumArea = cumsum(2*pi.*r_mid.*dose_r.*dr);
                      
                      if abs((idd(ixDepth(j)) * conversionFactor) - cumArea(end)) > 1e-3
                         warning('shell integration is wrong')
                      end
        
                      if cumArea(end) > iddPeak * (1-cutOffLevel)
                          IX = find(cumArea >= cumArea(end) * cutOffLevel,1, 'first'); 
                      else
                          IX = 1;
                      end
                      machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;
                case '3D'
                      warning('not yet implemented')
                      % find smallest enclosing volume that covers X-percent of the
                      % total integral dose
            end


            if isempty(IX)
                depthDoseCutOff =   0;
            elseif isnumeric(IX)
                depthDoseCutOff =  r_mid(IX);
            end

            machine.data(energyIx).LatCutOff.CutOff(j) = depthDoseCutOff;

        end
    end    
end    



            
%% visualization

if visBool
    
    % determine which pencil beam should be plotted
    subIx    = ceil(numel(vEnergiesIx)/2);
    energyIx = vEnergiesIx(subIx);
    
    baseData       = machine.data(energyIx);
    maxfocusIx     = energySigmaLUT(ix_Max(subIx),2);
    maxSSD         = energySigmaLUT(ix_Max(subIx),3);
    rangeShifter   = rangeShifterLUT(ix_Max(subIx));
    radiationMode  = stf(1).radiationMode;
    
    TmpCompFac     = baseData.LatCutOff.CompFac;
    baseData.LatCutOff.CompFac = 1;
    
    % plot 3D cutoff at one specific depth on a rather sparse grid
    sStep         = 0.25;
    vLatX         = -max(machine.data(energyIx).LatCutOff.CutOff)*3 : sStep : max(machine.data(energyIx).LatCutOff.CutOff)*3; % [mm]

    dimX          = numel(vLatX);
    midPos        = round(length(vLatX)/2);
    [X,Y]         = meshgrid(vLatX,vLatX);
    
    radDepths     = 0:1:machine.data(energyIx).depths(end);
    radialDist_sq = (X.^2 + Y.^2);
    radialDist_sq = radialDist_sq(:);
    mDose         = zeros(dimX,dimX,numel(radDepths));
    vDoseInt      = zeros(numel(radDepths),1);
    
    for kk = 1:numel(radDepths)          
         mDose(:,:,kk) = reshape(matRad_calcParticleDoseBixel(radDepths(kk)*ones(numel(radialDist_sq),1), radialDist_sq, maxSSD,...
              maxfocusIx, baseData, rangeShifter, radiationMode),[dimX dimX]);
          
         IX               =  find((radDepths(kk)-machine.data(energyIx).LatCutOff.depths)>=0,1,'last');
         TmpCutOff        = machine.data(energyIx).LatCutOff.CutOff(IX);
         
         vXCut = vX(vX<=TmpCutOff);
         
         % integration steps
         r_mid_Cut        = (0.5*(vXCut(1:end-1) +  vXCut(2:end)))'; % [mm]
         dr_Cut           = (vXCut(2:end) - vXCut(1:end-1))';
         radialDist_sqCut = r_mid_Cut.^2;
              
         dose_r_Cut    = matRad_calcParticleDoseBixel(radDepths(kk)*ones(numel(radialDist_sqCut),1), radialDist_sqCut(:), maxSSD,...
                                 maxfocusIx, baseData, rangeShifter, radiationMode);
         
         cumAreaCut = cumsum(2*pi.*r_mid_Cut.*dose_r_Cut.*dr_Cut);  
         
         if ~isempty(cumAreaCut)
             vDoseInt(kk) = cumAreaCut(end);
         end
    end
    
    % obtain maximum dose
    [~,peakixDepth] = max(machine.data(energyIx).Z); 
    dosePeakPos = matRad_calcParticleDoseBixel(machine.data(energyIx).depths(peakixDepth), 0, maxSSD, maxfocusIx, baseData, rangeShifter, radiationMode);   
    
    vLevelsDose = dosePeakPos.*[0.01 0.05 0.1 0.9];
    figure,set(gcf,'Color',[1 1 1]);
    subplot(121),h=imagesc(squeeze(mDose(midPos,:,:)));hold on;
    set(h,'AlphaData', .8*double(squeeze(mDose(fix(numel(vLatX)/2),:,:))>0));
    contour(squeeze(mDose(fix(numel(vLatX)/2),:,:)),vLevelsDose,'LevelListMode','manual','LineWidth',3);hold on
    plot(machine.data(energyIx).LatCutOff.depths, machine.data(energyIx).LatCutOff.CutOff * sStep^-1 + midPos,'rx');
    legend({'isodose 1%,5%,10% 90%','calculated cutoff'}) ,colorbar,set(gca,'FontSize',12),xlabel('z [mm]')',ylabel('x [mm]')
    subplot(122),plot(machine.data(energyIx).depths,machine.data(energyIx).Z*conversionFactor,'k','LineWidth',2),grid on,hold on
                 plot(radDepths,vDoseInt,'r--','LineWidth',2),hold on,
                 plot(radDepths,vDoseInt * TmpCompFac,'bx','LineWidth',1),hold on,
    legend({'original IDD',['cut off IDD at ' num2str(cutOffLevel) '%'],'cut off IDD with compensation'},'Location','northwest'),xlabel('z [mm]') ,set(gca,'FontSize',12)     
           
    totEnergy        = trapz(machine.data(energyIx).depths,machine.data(energyIx).Z*conversionFactor) ;
    totEnergyCutOff  = trapz(radDepths,vDoseInt * TmpCompFac) ;
    relDiff          =  ((totEnergy/totEnergyCutOff)-1)*100;   
    title(['rel diff of integral dose ' num2str(relDiff) '%']);

    baseData.LatCutOff.CompFac = TmpCompFac;
    
     % set depth position - 1 means plotting the entry profile
    j              = ixDepth(1);
    CutOff         = machine.data(energyIx).LatCutOff.CutOff(j);
   
    DoseSlice = mDose(:,:,j);
    [~,LevelixDepth] = min(abs(X(1,:)-ceil(CutOff)));
    DoseLevel = DoseSlice(midPos,LevelixDepth);
    
    figure,set(gcf,'Color',[1 1 1]);
    subplot(221),surf(X,Y,DoseSlice),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'), hold on, axis tight
    contour3(X,Y,DoseSlice,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r'),hold on;
    title({['beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)], ['cutoff = ' num2str(cutOffLevel)]}),set(gca,'FontSize',12);
    
    subplot(222),surf(X,Y,DoseSlice),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'), hold on, axis tight
    contour3(X,Y,DoseSlice,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r'),hold on; title(['intensity profile; cutoff = ' num2str(cutOffLevel)]),view(0,90)
    
    vDoseLat = mDose(round(numel(vLatX)/2),:,j);
    
    subplot(223),plot(vLatX,vDoseLat,'LineWidth',3),grid on, grid minor, hold on
    plot([CutOff,CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on
    plot([-CutOff,-CutOff],[0 max(vDoseLat)],'r','LineWidth',2),hold on, title('lateral profile 2D - cross section')
    
    subplot(224),surf(X,Y,DoseSlice),xlabel('x'),ylabel('y'),zlabel('double lateral gauss'),colormap(parula(256)), hold on
    title(['proton beam with energy ' num2str(machine.data(energyIx).energy) ' at depth index ' num2str(j)]),set(gca,'FontSize',12);
    contour3(X,Y,DoseSlice,[(DoseLevel+0.001*DoseLevel) DoseLevel],'LineWidth',3,'color','r');title('lateral profile 3D'); view([0 0]);


    % plot cutoff of different energies
    figure,set(gcf,'Color',[1 1 1]);
    cnt = 1;
    for i = vEnergiesIx
        plot(machine.data(i).LatCutOff.depths,machine.data(i).LatCutOff.CutOff,'LineWidth',1.5),hold on
        cellLegend{cnt} = [num2str(machine.data(i).energy) ' MeV'];
        cnt = cnt + 1;
    end
    grid on, grid minor,xlabel('depth in [mm]'),ylabel('lateral cutoff in [mm]')
    title(['cutoff level = ' num2str(cutOffLevel)]),
    ylim = get(gca,'Ylim');    set(gca,'Ylim',[0 ylim(2)+3]),    legend(cellLegend)
end





end

