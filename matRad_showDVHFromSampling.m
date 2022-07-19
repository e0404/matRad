function matRad_showDVHFromSampling(caSamp,cst,pln,scenarios,doseWindow,dvhType,refScen,lineStyleIndicator)
% matRad dvh visualizaion
%
% call
%   matRad_showDVH(dvh,cst,pln,lineStyleIndicator)
%
% input
%   result:             result struct from fluence optimization/sequencing
%   cst:                matRad cst struct
%   pln:                matRad pln struct
%   lineStyleIndicator: integer (1,2,3,4) to indicate the current linestyle
%                       (hint: use different lineStyles to overlay
%                       different dvhs)
%
% output
%   graphical display of DVH
%
% References
%   -
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

[env,envver] = matRad_getEnvironment();

if ~exist('dvhType','var') || isempty(dvhType)
    dvhType = 'multiscenario';
end

if ~exist('refScen','var') || isempty(refScen)
    refScen = 1;
end

if ~exist('lineStyleIndicator','var') || isempty(lineStyleIndicator)
    lineStyleIndicator = 1;
end

% create new figure and set default line style indicator if not explictly
% specified
hold on;

cstNames = cst(:,2);
cstInfo = cst(:,5);

numOfVois = numel(cstNames);


vProb = zeros(numel(caSamp),1);

for l = 1:numel(caSamp)
    
    [ctScen,shiftScen,RangeScen] = deal(pln.multScen.linearMask(l,1),pln.multScen.linearMask(l,2),pln.multScen.linearMask(l,3));
    shiftScenMask = find(squeeze(pln.multScen.scenMask(1,:,:)));
    indProb = sub2ind([pln.multScen.totNumShiftScen pln.multScen.totNumRangeScen],shiftScen,RangeScen);
    
    numCtScen = nnz(pln.multScen.scenMask(:,shiftScen,RangeScen));
    if(numCtScen>1)
        vProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb))/phaseProb(ctScen);
    else
        vProb(l)=pln.multScen.scenProb(find(shiftScenMask==indProb));
    end
end

%% print the dvh
try
    colorMx = cellfun(@(c) c.visibleColor,cstInfo,'UniformOutput',false);
    colorMx = cell2mat(colorMx);
catch
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/numOfVois):64,:);
end

lineStyles = {'-',':','--','-.'};

maxDVHvol  = 0;
maxDVHdose = 0;

switch dvhType
    
    case 'multiscenario'
        
        for i = 1:numOfVois
            for k = scenarios
                if cst{i,5}.Visible
                    % cut off at the first zero value where there is no more signal
                    % behind
                    ix      = max([1 find(caSamp(k).dvh(i).volumePoints>0,1,'last')]);
                    currDvh = [caSamp(k).dvh(i).doseGrid(1:ix);caSamp(k).dvh(i).volumePoints(1:ix)];
                    
                    p=plot(currDvh(1,:),currDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                        'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});
                    
                    if(k==refScen)
                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    else
                        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end
                    
                    maxDVHvol  = max(maxDVHvol,max(currDvh(2,:)));
                    maxDVHdose = max(maxDVHdose,max(currDvh(1,:)));
                    
                end
            end
        end
        
    case 'minmax'
        
        for i = 1:numOfVois
            
            allDvh = caSamp(1).dvh(1).doseGrid;
            
            for k = scenarios
                if cst{i,5}.Visible
                    if(k==refScen)
                        % cut off at the first zero value where there is no more signal
                        % behind
                        ix      = max([1 find(caSamp(k).dvh(i).volumePoints>0,1,'last')]);
                        refDvh = [caSamp(k).dvh(i).doseGrid(1:ix);caSamp(k).dvh(i).volumePoints(1:ix)];
                        
                        p=plot(refDvh(1,:),refDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});
                        
                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    end
                    
                    allDvh = [allDvh;caSamp(k).dvh(i).volumePoints];
                end
            end
            
            bandDvh = [allDvh(1,:);min(allDvh(2:end,:));max(allDvh(2:end,:))];
            
            f = fill([bandDvh(1,:) flip(bandDvh(1,:))],[bandDvh(2,:) flip(bandDvh(3,:))],colorMx(i,:));
            f.FaceAlpha = 0.2;
            f.FaceColor = colorMx(i,:);
            f.EdgeColor = colorMx(i,:);
            f.LineWidth = 0.05;
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            maxDVHvol  = max(maxDVHvol,max(bandDvh(3,:)));
            maxDVHdose = max(maxDVHdose,max(bandDvh(1,:)));
            
        end
        
    case 'trustband'
        
        for i = 1:numOfVois
            
            allDvh = caSamp(1).dvh(1).doseGrid;
            
            for k = scenarios
                if cst{i,5}.Visible
                    if(k==refScen)
                        % cut off at the first zero value where there is no more signal
                        % behind
                        ix      = max([1 find(caSamp(k).dvh(i).volumePoints>0,1,'last')]);
                        refDvh = [caSamp(k).dvh(i).doseGrid(1:ix);caSamp(k).dvh(i).volumePoints(1:ix)];
                        
                        p=plot(refDvh(1,:),refDvh(2,:),'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                            'LineStyle',lineStyles{lineStyleIndicator},'DisplayName',cst{i,2});
                        
                        p.LineWidth = 2;
                        p.Color = [colorMx(i,:) 1];
                        p.Annotation.LegendInformation.IconDisplayStyle = 'on';
                    end
                    
                    allDvh = [allDvh;caSamp(k).dvh(i).volumePoints];
                end
            end
                       
            meanDVHVol  = mean(allDvh(2:end,:),1);
            stdDVHVol   = std(allDvh(2:end,:),1,1);
            meanDVHVolW = (sum(allDvh(2:end,:)' * diag(vProb),2))';
            
            %Weighting of std is not similar in Octave & Matlab
            switch env
                case 'MATLAB'              
                    stdDVHVolW  = std(allDvh(2:end,:),vProb,1);
                    
                case 'OCTAVE'
                    try
                        pkg load nan;
                        nanLoaded = true;
                    catch
                        nanLoaded = false;
                        matRad_cfg.dispWarning('Weighted std not possible due to missing ''nan'' package!');
                    end
                    
                    if nanLoaded
                        stdDVHVolW  = std(allDvh(2:end,:),[],1,vProb);
                    else
                        stdDVHVolW  = stdCube;
                    end
            end
            
            bandDvh = [allDvh(1,:);meanDVHVolW-stdDVHVolW;meanDVHVolW+stdDVHVolW];
            
            p=plot(bandDvh(1,:),meanDVHVolW,'LineWidth',0.5,'Color',[colorMx(i,:) 0.4], ...
                'LineStyle',lineStyles{2},'DisplayName',cst{i,2});
            
            p.LineWidth = 2;
            p.Color = [colorMx(i,:) 1];
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            f = fill([bandDvh(1,:) flip(bandDvh(1,:))],[bandDvh(2,:) flip(bandDvh(3,:))],colorMx(i,:));
            f.FaceAlpha = 0.2;
            f.FaceColor = colorMx(i,:);
            f.EdgeColor = colorMx(i,:);
            f.LineWidth = 0.05;
            f.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            maxDVHvol  = max(maxDVHvol,max([max(bandDvh(2,:)) max(bandDvh(3,:))]));
            maxDVHdose = max(maxDVHdose,max(bandDvh(1,:)));
            
        end
        
end

fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',10,'Interpreter','none');
legend boxoff

if ~exist('doseWindow', 'var')
    doseWindow = [0 1.4*maxDVHdose];
end

ylim([0 1.05*maxDVHvol]);
xlim(doseWindow);

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',.5,'FontSize',fontSizeValue);
ylabel('Volume [%]','FontSize',fontSizeValue)

if strcmp(pln.bioParam.model,'none')
    xlabel('Dose [Gy]','FontSize',fontSizeValue);
else
    xlabel('RBE x Dose [Gy(RBE)]','FontSize',fontSizeValue);
end