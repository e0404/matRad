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

if cutOffLevel <= 0.98
   warning('a lateral cut off below 0.98 may result in an inaccurate dose calculation') 
end

TypeOfCutOffCalc = '2D'; % 'classic','2D','3D',
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
vX = [0 logspace(-1,3,1200)]; % [mm]

% integration steps
r_mid          = (0.5*(vX(1:end-1) +  vX(2:end)))'; % [mm]
dr             = (vX(2:end) - vX(1:end-1))';
radialDist_sq  = r_mid.^2;

% number of depth points for which a lateral cutoff is determined
NumDepthVal       = 35; 
epsilon           = 1e-5;  

% define function handles for single and double gauss
CF    =  @(LcutSigma)(1/(1-exp((-LcutSigma^2)/2)));
% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

% extract SSD for each bixel
vSSD = ones(1,length([stf.ray(:).energy]));
cnt = 1;
for i  = 1:length(stf.ray)
    vSSD(cnt:cnt+numel([stf.ray(i).energy])-1) = stf.ray(i).SSD;
    cnt = cnt + numel(stf.ray(i).energy);
end

% setup energy, focus index, sigma look up table - only consider unique rows
[energySigmaLUT,ixUnique]  = unique([[stf.ray(:).energy]; [stf.ray(:).focusIx] ; vSSD]','rows');
rangeShifterLUT = [stf.ray(:).rangeShifter];
rangeShifterLUT = rangeShifterLUT(1,ixUnique);

% find the largest inital beam width considering focus index, SSD and range shifter for each individual energy
for i = 1:size(energySigmaLUT,1)
   
    % find index of maximum used energy (round to keV for numerical reasons
    energyIx = max(round2(energySigmaLUT(i,1),4)) == round2([machine.data.energy],4);
    
    currFoci = energySigmaLUT(i,2);
    sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist(currFoci,:)',...
                              machine.data(energyIx).initFocus.sigma(currFoci,:)',...
                              energySigmaLUT(i,3));
    sigmaIni_sq = sigmaIni^2;
    
    % consider range shifter for protons if applicable
    if  strcmp(machine.meta.radiationMode,'protons') && rangeShifterLUT(i).eqThickness > 0  && ~strcmp(machine.meta.machine,'Generic')

        %get max range shift
        sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy, ...
                                           rangeShifterLUT(i), ...
                                           energySigmaLUT(i,3));

        % add to initial sigma in quadrature
        sigmaIni_sq = sigmaIni_sq +  sigmaRashi.^2;

    end                          
                                                         
    energySigmaLUT(i,4) = sigmaIni_sq;
    
end

% find for each individual energy the broadest inital beam width
uniqueEnergies                = unique(energySigmaLUT(:,1));
largestSigmaSq4uniqueEnergies = NaN * ones(numel(uniqueEnergies),1);
ix_Max                        = NaN * ones(numel(uniqueEnergies),1);
for i = 1:numel(uniqueEnergies)
    [largestSigmaSq4uniqueEnergies(i), ix_Max(i)] = max(energySigmaLUT(uniqueEnergies(i) == energySigmaLUT(:,1),4));
end

% get energy indices for looping
vEnergiesIx = find(ismember([machine.data(:).energy],uniqueEnergies(:,1)));
cnt         = 0;    

% loop over all entries in the machine.data struct
for energyIx = vEnergiesIx
   
    % set default depth cut off - finite value will be set during first iteration
    depthDoseCutOff = inf;

    % get the current integrated depth dose profile
    if isstruct(machine.data(energyIx).Z)
        idd_org = SumGauss(machine.data(energyIx).depths,machine.data(energyIx).Z.mean,...
                                   machine.data(energyIx).Z.width.^2,...
                                   machine.data(energyIx).Z.weight) * conversionFactor;
    else
        idd_org = matRad_interp1(machine.data(energyIx).depths,machine.data(energyIx).Z,machine.data(energyIx).depths) * conversionFactor;
    end
    
    [~,peakIxOrg] = max(idd_org); 
    
    % get indices for which a lateral cutoff should be calculated
    cumIntEnergy = cumtrapz(machine.data(energyIx).depths,idd_org);
    
    if strcmp(machine.meta.radiationMode,'protons')
        vEnergySteps        = 0:(cumIntEnergy(end)/NumDepthVal):cumIntEnergy(end);
    else
        peakTailRelation   = 0.5;
        NumDepthValToPeak  = ceil(NumDepthVal*peakTailRelation);                                                                          % number of depth values from 0 to peak position
        NumDepthValTail    = ceil(NumDepthVal*(1-peakTailRelation));                                                                      % number of depth values behind peak position
        EnergyStepsToPeak  = cumIntEnergy(peakIxOrg)/NumDepthValToPeak;
        EnergyStepsTail    = (cumIntEnergy(end)-cumIntEnergy(peakIxOrg))/NumDepthValTail;
        vEnergySteps       = [0:EnergyStepsToPeak:cumIntEnergy(peakIxOrg) cumIntEnergy(peakIxOrg+1):EnergyStepsTail:cumIntEnergy(end)];
    end
    
    [cumIntEnergy,ix] = unique(cumIntEnergy);
    depthValues       = matRad_interp1(cumIntEnergy,machine.data(energyIx).depths(ix),vEnergySteps);
    idd               = matRad_interp1(machine.data(energyIx).depths,idd_org,depthValues);          
    [~,peakIx]        = max(idd); 
    
    cnt = cnt +1 ;
    % % calculate dose in spot
    baseData                   = machine.data(energyIx);
    baseData.LatCutOff.CompFac = 1;   
  
    for j = 1:numel(depthValues)
        
        % save depth value
        machine.data(energyIx).LatCutOff.depths(j) = depthValues(j);
        
        radDepths      = (depthValues(j) + baseData.offset - epsilon) * ones(numel(r_mid),1); 
        
        dose_r         = matRad_calcParticleDoseBixel(depthValues(j) + baseData.offset, radialDist_sq, largestSigmaSq4uniqueEnergies(cnt), baseData);
       
             
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
    
                case '2D'
                      cumArea = cumsum(2*pi.*r_mid.*dose_r.*dr);
                      relativeThreshold = 0.5; %in [%]
                      if abs((cumArea(end)./(idd(j)))-1)*100 > relativeThreshold && (cumArea(end) > relativeThreshold * idd(peakIx))
                         warning('LateralParticleCutOff: shell integration is wrong')
                      end
        
                      IX = find(cumArea >= idd(j) * cutOffLevel,1, 'first'); 
                      machine.data(energyIx).LatCutOff.CompFac = cutOffLevel^-1;
                      
                case '3D'
                      warning('LateralParticleCutOff: not yet implemented')
                      % find smallest enclosing volume that covers X-percent of the total integral dose
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
    focusIx        = energySigmaLUT(ix_Max(subIx),2);
    maxSSD         = energySigmaLUT(ix_Max(subIx),3);
    rangeShifter   = rangeShifterLUT(ix_Max(subIx));
    TmpCompFac     = baseData.LatCutOff.CompFac;
    baseData.LatCutOff.CompFac = 1;
    
    % plot 3D cutoff at one specific depth on a rather sparse grid
    sStep         = 0.5;
    vLatX         = -100 : sStep : 100; % [mm]
    dimX          = numel(vLatX);
    midPos        = round(length(vLatX)/2);
    [X,Y]         = meshgrid(vLatX,vLatX);
    
    radDepths     = [0:sStep:machine.data(energyIx).depths(end)] + machine.data(energyIx).offset;
    radialDist_sq = (X.^2 + Y.^2);
    radialDist_sq = radialDist_sq(:);
    mDose         = zeros(dimX,dimX,numel(radDepths));
    vDoseInt      = zeros(numel(radDepths),1);
    
    for kk = 1:numel(radDepths)    
          
         % calculate initial focus sigma
         sigmaIni = matRad_interp1(machine.data(energyIx).initFocus.dist(focusIx,:)', ...
                                   machine.data(energyIx).initFocus.sigma(focusIx,:)',maxSSD);
         sigmaIni_sq = sigmaIni^2;

         % consider range shifter for protons if applicable
         if rangeShifter.eqThickness > 0 && strcmp(pln.radiationMode,'protons')

              % compute!
              sigmaRashi = matRad_calcSigmaRashi(machine.data(energyIx).energy,rangeShifter,maxSSD);

              % add to initial sigma in quadrature
              sigmaIni_sq = sigmaIni_sq +  sigmaRashi^2;

         end
                      
       
         mDose(:,:,kk) = reshape(matRad_calcParticleDoseBixel(radDepths(kk), radialDist_sq, sigmaIni_sq,baseData),[dimX dimX]);
          
         [~,IX]           = min(abs((machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset) - radDepths(kk)));
         TmpCutOff        = machine.data(energyIx).LatCutOff.CutOff(IX);    
         vXCut            = vX(vX<=TmpCutOff);
         
         % integration steps
         r_mid_Cut        = (0.5*(vXCut(1:end-1) +  vXCut(2:end)))'; % [mm]
         dr_Cut           = (vXCut(2:end) - vXCut(1:end-1))';
         radialDist_sqCut = r_mid_Cut.^2;    
         
         dose_r_Cut       = matRad_calcParticleDoseBixel(radDepths(kk), radialDist_sqCut(:), sigmaIni_sq,baseData);
         
         cumAreaCut = cumsum(2*pi.*r_mid_Cut.*dose_r_Cut.*dr_Cut);  
         
         if ~isempty(cumAreaCut)
             vDoseInt(kk) = cumAreaCut(end);
         end
    end
    
    % obtain maximum dose
    [~,peakixDepth] = max(machine.data(energyIx).Z); 
    dosePeakPos = matRad_calcParticleDoseBixel(machine.data(energyIx).depths(peakixDepth), 0, sigmaIni_sq, baseData);   
    
    vLevelsDose = dosePeakPos.*[0.01 0.05 0.1 0.9];
    doseSlice   = squeeze(mDose(midPos,:,:));
    figure,set(gcf,'Color',[1 1 1]);
    subplot(311),h=imagesc(squeeze(mDose(midPos,:,:)));hold on;
    set(h,'AlphaData', .8*double(doseSlice>0));
    contour(doseSlice,vLevelsDose,'LevelListMode','manual','LineWidth',2);hold on
    
    ax = gca;
    ax.XTickLabelMode = 'manual';
    ax.XTickLabel     = strsplit(num2str(ax.XTick*sStep + machine.data(energyIx).offset),' ')';
    ax.YTickLabelMode = 'manual';
    ax.YTickLabel     = strsplit(num2str(ax.YTick*sStep + machine.data(energyIx).offset),' ')';
  
    plot(1+(machine.data(energyIx).LatCutOff.depths)*sStep^-1,...
          machine.data(energyIx).LatCutOff.CutOff * sStep^-1 + midPos,'rx');

    legend({'isodose 1%,5%,10% 90%','calculated cutoff'}) ,colorbar,set(gca,'FontSize',12),xlabel('z [mm]'),ylabel('x [mm]');
       
    entry = machine.data(energyIx);
    if isstruct(entry.Z)
       idd = SumGauss(entry.depths,entry.Z.mean,entry.Z.width.^2,entry.Z.weight);
    else
        idd = machine.data(energyIx).Z;
    end
    subplot(312),plot(machine.data(energyIx).depths,idd*conversionFactor,'k','LineWidth',2),grid on,hold on
                 plot(radDepths - machine.data(energyIx).offset,vDoseInt,'r--','LineWidth',2),hold on,
                 plot(radDepths - machine.data(energyIx).offset,vDoseInt * TmpCompFac,'bx','LineWidth',1),hold on,
    legend({'original IDD',['cut off IDD at ' num2str(cutOffLevel) '%'],'cut off IDD with compensation'},'Location','northwest'),
    xlabel('z [mm]'),ylabel('[MeV cm^2 /(g * primary)]'),set(gca,'FontSize',12)     
           
    totEnergy        = trapz(machine.data(energyIx).depths,machine.data(energyIx).Z*conversionFactor) ;
    totEnergyCutOff  = trapz(radDepths,vDoseInt * TmpCompFac) ;
    relDiff          =  ((totEnergy/totEnergyCutOff)-1)*100;   
    title(['rel diff of integral dose ' num2str(relDiff) '%']);
    baseData.LatCutOff.CompFac = TmpCompFac;
    
    subplot(313),
    if isfield(machine.data(energyIx),'sigma1')
        yyaxis left;
        plot(machine.data(energyIx).LatCutOff.depths,machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on
        plot(machine.data(energyIx).depths,(machine.data(energyIx).sigma1),':','LineWidth',2),grid on,hold on,ylabel('mm')
        yyaxis right; 
        plot(machine.data(energyIx).depths,(machine.data(energyIx).sigma2),'-.','LineWidth',2),grid on,hold on,ylabel('mm')
        legend({'Cutoff','sigma1','sigma2'});
    else
        yyaxis left;plot(machine.data(energyIx).LatCutOff.depths,machine.data(energyIx).LatCutOff.CutOff,'LineWidth',2),hold on,ylabel('mm')
        yyaxis right;subplot(313),plot(machine.data(energyIx).depths,machine.data(energyIx).sigma,'LineWidth',2),grid on,hold on
        legend({'Cutoff','sigma'});ylabel('mm')
    end

    set(gca,'FontSize',12),xlabel('z [mm]'),  ylabel('mm')

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

