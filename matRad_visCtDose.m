 function matRad_visCtDose(optResult,cst,pln,ct,slice,plane)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad visualization of two-dimensional dose distributions on ct including
% segmentation
% 
% call
%   matRad_visCtDose(dose,cst,pln,ct,slice,plane)
%
% input
%   dose:   dose cube
%   cst:    matRad cst struct
%   pln:    matRad plan meta information struct
%   ct:     ct cube
%   slice:  number of slice for initial visualization (optional)
%   plane:  integer to select coronal (1), sagital (2), or axial (3) plane (optional)
%
% output
%   -   
%
% References
%   -
%
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

if nargin > 0
    data.optResult = optResult;
    data.cst  = cst;
    data.pln  = pln;
    data.ct   = ct;
    if ~isempty(data.optResult)
        data.optResult = rmfield(data.optResult,'w');
        data.fName =fieldnames(data.optResult);
        for i=1:size(data.fName,1)
            %indicates if it should be plotted later on
            data.fName{i,2}=1;
            % Reshape dose to cube in voxel dimensions
            CurrentCube = getfield(data.optResult,data.fName{i,1});
            if ~isempty(CurrentCube) && ~isempty(data.ct) && isequal(size(CurrentCube),size(data.ct))
                data.optResult = setfield(data.optResult,data.fName{i,1},reshape(CurrentCube,size(ct)));
            %try to reshape using voxelDimensions from pln struct    
            elseif ~isempty(data.optResult) && ~isempty(data.pln) && isequal(size(CurrentCube),size(data.ct))
                data.optResult = setfield(data.optResult,data.fName{i,1},reshape(CurrentCube,data.pln.voxelDimensions));
            elseif ~isempty(data.optResult) && ~strcmp(data.fName{i,1},'w')
                error('Cannot reshape dose');   
            end
        end
    end
    
    if ~isempty(data.optResult)
        data.doseColorwashCheckboxValue = 1;
        data.doseIsoCheckboxValue = 1;
        data.SelectedDisplayOption = 'physicalDose';
        data.TypeOfPlot = 1;
    else
        data.doseColorwashCheckboxValue = 0;
        data.doseIsoCheckboxValue = 0;
        data.SelectedDisplayOption = 'physicalDose';
        data.TypeOfPlot = 1;
    end
    
    if ~isempty(data.ct)
        data.ctCheckboxValue = 1;
    else
        data.ctCheckboxValue = 0;
    end
    
    if ~isempty(data.cst)
        data.contourCheckboxValue = 1;
    else
        data.contourCheckboxValue = 0;
    end
        
        % set standard values
    if nargin < 6
        data.plane = 3;
    else
        data.plane = plane;
    end
    
    if nargin < 5
        data.slice = round(pln.isoCenter(data.plane)/pln.resolution(data.plane));
    else
        data.slice = slice;
    end
    
    if data.plane == 1
        data.axis = [1 size(data.ct,1) 1 size(data.ct,3)];
    elseif data.plane == 2        
        data.axis = [1 size(data.ct,3) 1 size(data.ct,2)];
    elseif data.plane == 3
        data.axis = [1 size(data.ct,1) 1 size(data.ct,2)];
    end
    
    data.LateralOffset = NaN;
    
    % Open figure
    myWindow = figure('Name','matRad CT/dose/VOI bowser','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'ToolBar','figure');
    myAxes   = axes('Position', [0.35 0.1 0.55 0.8],'YDir','reverse');
    
    guidata(myWindow,data);
    
else
    
    myWindow = gcf;
    data = guidata(myWindow);
    
    vAxis = [1 size(data.ct,1) 1 size(data.ct,2)];
    data.axis = vAxis;
    
    clf;
    
    myAxes   = axes('Position', [0.35 0.1 0.55 0.8],'YDir','reverse');
        
end

    hold on


    %% ct
if ~isempty(data.ct) && data.ctCheckboxValue && data.TypeOfPlot ==1
    
    if data.plane == 1 % Coronal plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct(data.slice,:,:))/max(data.ct(:))),bone);
    elseif data.plane == 2 % Sagital plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct(:,data.slice,:))/max(data.ct(:))),bone);
    elseif data.plane == 3 % Axial plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct(:,:,data.slice))/max(data.ct(:))),bone);
    end
    
    axes(myAxes)
    ctImageHandle = image(ct_rgb);

end

if ~isempty(data.optResult) && data.TypeOfPlot ==1
    mVolume = getfield(data.optResult,data.SelectedDisplayOption);
%     %% dose colorwash
    if ~isempty(mVolume) && data.doseColorwashCheckboxValue && ~isvector(mVolume)

        dose_rgb = mVolume./max(mVolume(:));

        % Save RGB indices for dose in zsliceÂ´s voxels.
        if data.plane == 1  % Axial plane
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(data.slice,:,:))),jet);
        elseif data.plane == 2 % Sagital plane
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,data.slice,:))),jet);
        elseif data.plane == 3 % Coronal plane
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,data.slice))),jet);
        end
        % Show dose
        axes(myAxes);
        doseImageHandle = image(dose_rgb);

        % Make dose transparent
        if ~isempty(data.ct)
            set(doseImageHandle,'AlphaData',.45);
        end

    end

    %% dose iso dose lines
    if ~isempty(mVolume) && data.doseIsoCheckboxValue && data.TypeOfPlot ==1 && ~isvector(mVolume)
         
        delta = (max(mVolume(:))-min(mVolume(:)))*0.1;
        vSpacingIsoDose = linspace(min(mVolume(:))+delta,max(mVolume(:))+delta,10);
        vSpacingIsoDose= round(vSpacingIsoDose.*100)/100;
       % vSpacingIsoDose = 5:5:100;
        
        if data.plane == 1  % Coronal plane
            [~,myContour] = contour(myAxes,squeeze(mVolume(data.slice,:,:)),vSpacingIsoDose);
        elseif data.plane == 2 % Sagittal plane
            [~,myContour] = contour(myAxes,squeeze(mVolume(:,data.slice,:)),vSpacingIsoDose);
        elseif data.plane == 3 % Axial plane
             [~,myContour] = contour(myAxes,squeeze(mVolume(:,:,data.slice)),vSpacingIsoDose);
        end

        % turn off legend for this data set
        hAnnotation = get(myContour,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')

        set(myContour,'LabelSpacing',80,'ShowText','on')
    end

end
    %% VOIs
if ~isempty(data.cst) && data.contourCheckboxValue && data.TypeOfPlot ==1

    colors = jet;
    colors = colors(round(linspace(1,63,size(data.cst,1))),:);

    mask = zeros(size(data.ct)); % create zero cube with same dimeonsions like dose cube
    for s = 1:size(data.cst,1)
        if ~strcmp(data.cst{s,3},'IGNORED') %&& ~strcmp(data.cst{s,2},'DoseFalloff')
            mask(:) = 0;
            mask(data.cst{s,8}) = 1;
            if data.plane == 1 && sum(sum(mask(data.slice,:,:))) > 0
                contour(myAxes,squeeze(mask(data.slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            elseif data.plane == 2 && sum(sum(mask(:,data.slice,:))) > 0
                contour(myAxes,squeeze(mask(:,data.slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            elseif data.plane == 3 && sum(sum(mask(:,:,data.slice))) > 0
                contour(myAxes,squeeze(mask(:,:,data.slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            end
        end
    end

    myLegend = legend('show');
    set(myLegend,'Position',[.85 .5 .1 .2]);
    set(myLegend,'FontSize',20);
    
end
    
%% Set axis labels
if   data.plane == 3% Axial plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.pln.resolution(1):1000)
        set(gca,'YTick',0:50/data.pln.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*data.pln.resolution(1))
        set(gca,'YTickLabel',0:50:1000*data.pln.resolution(2))
        xlabel('x [mm]')
        ylabel('y [mm]')
        title('Axial plane')
    else
        xlabel('x [voxels]')
        ylabel('y [voxels]')
        title('Axial plane')
    end
elseif data.plane == 2 % Sagittal plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.pln.resolution(3):1000)
        set(gca,'YTick',0:50/data.pln.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*data.pln.resolution(3))
        set(gca,'YTickLabel',0:50:1000*data.pln.resolution(2))
        xlabel('z [mm]')
        ylabel('y [mm]')
        title('Sagital plane');
    else
        xlabel('z [voxels]')
        ylabel('y [voxels]')
        title('Sagital plane');
    end
elseif data.plane == 1 % Coronal plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.pln.resolution(3):1000)
        set(gca,'YTick',0:50/data.pln.resolution(1):1000)
        set(gca,'XTickLabel',0:50:1000*data.pln.resolution(3))
        set(gca,'YTickLabel',0:50:1000*data.pln.resolution(1))
        xlabel('z [mm]')
        ylabel('x [mm]')
        title('Coronal plane')
    else
        xlabel('z [voxels]')
        ylabel('x [voxels]')
        title('Coronal plane')
    end
end

axis equal;
axis(data.axis);
    


if data.TypeOfPlot ==2

    % clear view and initialize some values
    clf;
    myAxes = axes('Position', [0.35 0.1 0.55 0.8]);
    set(gca,'YDir','normal');
    ylabel('dose')
    cColor={'k','g','m','y','c','r','b'};
    ymax=0;
    %a certain width to the profile plot can be added
    delta =0; % make it bixel distance dependend
    Cnt=1;
    
    mPhysDose=getfield(data.optResult,'physicalDose');
    mRotActualSlice =imrotate(mPhysDose(:,:,data.slice),data.pln.gantryAngles(1),'crop');
    
    
    vW =ones(size(mRotActualSlice,2),1);
    vProjected =vW'*mRotActualSlice;
    %find first and last nonzero element of projected data 
    %which can be seen as determinng the beam width
    idx = find(vProjected);
    % calculated average index to asses central axis
    idxCentAxis = round((idx(end)+idx(1))/2);
    % use central axis index or index from slider
    if isnan(data.LateralOffset)
        data.LateralOffset = idxCentAxis;
    else
        idxCentAxis=data.LateralOffset;
    end
    
    % plot physical dose
    mY=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
    mY(isnan(mY))=0;
    mY_avg=mean(mY,2);
    vX=linspace(1,data.pln.resolution(1)*numel(mY_avg),numel(mY_avg));
    PlotHandles{1} = plot(vX,mY_avg,'color',cColor{1,1},'LineWidth',3); hold on; 
    PlotHandles{1,2}='physicalDose';
    set(gca,'FontSize',18);
    % assess x - limits
    xLim  = find(mY_avg);
    if ~isempty(xLim)
        xmin= xLim(1)*data.pln.resolution(1)-20;
        xmax= xLim(end)*data.pln.resolution(1)+20;
    else
        vLim = axis;
        xmin = vLim(1);
        xmax = vLim(2);
    end
    
    % plot counter
    Cnt=2;
    
    if data.pln.bioOptimization == 1
        
        %will disable beta-plot
        %data.fName{7,2}=0;
        %will disable alpha-plot
        %data.fName{6,2}=0;
        
        for i=1:1:length(data.fName)
            mCurrentCube = getfield(data.optResult,data.fName{i,1});
            if ~isvector(mCurrentCube) && ~strcmp(data.fName{i,1},'RBEWeightedDose') ...
                    && ~strcmp(data.fName{i,1},'RBE') && ~strcmp(data.fName{i,1},'physicalDose')...
                    && data.fName{i,2}
                mRotActualSlice = imrotate(mCurrentCube(:,:,data.slice),data.pln.gantryAngles(1),'crop');
                mY = mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
                mY(isnan(mY))=0;
                mY=mean(mY,2);
                PlotHandles{Cnt,1} = plot(vX,mY,'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                PlotHandles{Cnt,2} =data.fName{i,1};
                Cnt = Cnt+1;
            end           
        end
        
        % plot always RBEWeightedDose against RBE
        mRBEWeightedDose=getfield(data.optResult,'RBEWeightedDose');
        mRotActualSlice =imrotate(mRBEWeightedDose(:,:,data.slice),data.pln.gantryAngles(1),'crop');
        mBED=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mBED(isnan(mBED))=0;
        vBED=mean(mBED,2);
        
        mRBE=getfield(data.optResult,'RBE');
        mRotActualSlice =imrotate(mRBE(:,:,data.slice),data.pln.gantryAngles(1),'crop');
        mRBE=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mRBE(isnan(mRBE))=0;
        vRBE=mean(mRBE,2);
        
        % plot biological dose against RBE
        [ax, PlotHandles{Cnt,1}, PlotHandles{Cnt+1,1}]=plotyy(vX,vBED,vX,vRBE,'plot');hold on;
        PlotHandles{Cnt,2}='RBEWeightedDose';
        PlotHandles{Cnt+1,2}='RBE';
         
        % set plotyy properties
        set(get(ax(2),'Ylabel'),'String','RBE','FontSize',18);
        set(get(ax(1),'Ylabel'),'String','RBE x dose','FontSize',18);
        set(PlotHandles{Cnt,1},'Linewidth',4,'color','r');
        set(PlotHandles{Cnt+1,1},'Linewidth',3,'color','b');
        set(ax(1),'ycolor','r')
        set(ax(2),'ycolor','b')
        set(ax,'FontSize',18);
        Cnt=Cnt+1;
       
    end
       
    
    % asses the prescripted dose and target coordinates 
    % todo: ptv, ctv, gtv are all labeld as target -> priorities
    sPrescrpDose=0;
    for i=1:size(data.cst,1)
        if strcmp(data.cst{i,3},'TARGET')==1
           mTarget = unique(data.cst{i,8});
           sPrescrpDose = data.cst{i,4};
        end
    end
    
     % plot prescription
    if sum(strcmp(fieldnames(data.optResult),'RBEWeightedDose')) > 0
            sPrescrpDose = sPrescrpDose./data.pln.numFractions;
    end
    PlotHandles{Cnt,1}=plot([0 size(data.ct,1)*data.pln.resolution(1)],[sPrescrpDose sPrescrpDose],'--','Linewidth',3,'color','m');
    PlotHandles{Cnt,2}='prescription';
    str = sprintf('profile plot of zentral axis of first beam at %d° at %d / %d in slice %d',data.pln.gantryAngles(1),data.LateralOffset*data.pln.resolution(2),size(data.ct,2)*data.pln.resolution(2), data.slice);
    title(str,'FontSize',14),grid on
    Cnt = Cnt+1;
    
    
    
    % plot target boundaries
    mTargetStack = zeros(size(data.ct));
    mTargetStack(mTarget)=1;
    mRotTargetSlice =imrotate(mTargetStack(:,:,data.slice),data.pln.gantryAngles(1),'crop');
    vRay = find(mRotTargetSlice(:,idxCentAxis))*data.pln.resolution(2);
    
    PlotHandles{Cnt,2} ='target boundary';
    vLim = axis;
    if ~isempty(vRay)
        PlotHandles{Cnt,1}=plot([vRay(1) vRay(1)],[0 vLim(4)],'--','Linewidth',2,'color','k');hold on
        plot([vRay(end) vRay(end)], [0 vLim(4)],'--','Linewidth',2,'color','k');hold on
        xmax = vRay(end)+30;
    else
        PlotHandles{Cnt,1} =0;
    end
    

    legend([PlotHandles{:,1}],PlotHandles{:,2},'Location','NorthWest');
    
    
    % set axis limits
    if data.pln.bioOptimization == 0 || ~strcmp(data.pln.radiationMode,'carbon')
        xlim([xmin xmax]);
    else
        xlim(ax(1),[xmin xmax]);
        xlim(ax(2),[xmin xmax]);
    end
    xlabel('depth [mm]','FontSize',16);
   
    

end
    
%% definition of ui
dataSetText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Data sets',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.04 0.86 0.06 0.03]);

doseColorwashCheckbox = uicontrol('Parent', gcf,...
        'Style', 'checkbox',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Dose (Colorwash)',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.05 0.82 0.096 0.03],...
        'Callback', @doseColorwashCheckboxCallback,...
        'Value',data.doseColorwashCheckboxValue);

doseIsoCheckbox = uicontrol('Parent', gcf,...
        'Style', 'checkbox',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Dose (Isodose lines)',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.05 0.78 0.096 0.03],...
        'Callback', @doseIsoCheckboxCallback,...
        'Value',data.doseIsoCheckboxValue);
           
ctCheckbox = uicontrol('Parent', gcf,...
        'Style', 'checkbox',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'CT',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.05 0.74 0.096 0.03],...
        'Callback', @ctCheckboxCallback,...
        'Value',data.ctCheckboxValue);
    
contourCheckbox = uicontrol('Parent', gcf,...
        'Style', 'checkbox',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Contour',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.05 0.7 0.096 0.03],...
        'Callback', @contourCheckboxCallback,...
        'Value',data.contourCheckboxValue);

    dataSetText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Plane',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.03 0.24 0.06 0.03]);

planePopup = uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', {'coronal','sagittal','axial'},...
        'FontSize', 10,...
        'Value',data.plane,...
        'Units', 'normalized',...
        'Position', [0.05 0.2 0.109 0.03],...
        'Callback', @planepopupCallback);

dataSetText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String',['Slice ' num2str(data.slice)],...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.03 0.14 0.06 0.03]);
    
slider = uicontrol('Parent', gcf,...
        'Style', 'slider',...
        'Units', 'normalized',...
        'Position', [0.05 0.10 0.2 0.02],...
        'Min', 1,...
        'Max',  size(data.ct,data.plane),...
        'Value', data.slice,...
        'SliderStep',[1/(size(data.ct,data.plane)-1) 1/(size(data.ct,data.plane)-1)],...
        'Callback', @sliderCallback);
    
doseSetText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Display options:',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.03 0.43 0.109 0.03],...
        'FontSize',14);    ;
 
if isempty(data.optResult)
    strTmp ={'no options available'};
else
    strTmp = data.fName(:,1);
    idx = find(strcmp(data.fName(:,1),data.SelectedDisplayOption));
end

dosePopup = uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', strTmp ,...
        'FontSize', 10,...
        'Value',idx,...
        'Units', 'normalized',...
        'Position', [0.05 0.4 0.109 0.03],...
        'Callback', @dosepopupCallback);

PopUpTypeOfPlotText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Type of plot:',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.025 0.53 0.109 0.03],...
        'FontSize',14);    
    
PopUpTypeOfPlot = uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', {'Intensity Plot','Profile Plot'} ,...
        'FontSize', 10,...
        'Value',data.TypeOfPlot,...
        'Units', 'normalized',...
        'Position', [0.05 0.50 0.059 0.03],...
        'Callback', @displaypopupCallback);

if data.TypeOfPlot == 2
       TmpStr ='on';
else
       TmpStr ='off';
end
ProfileSlider = uicontrol('Parent', gcf,...
        'Style', 'slider',...
        'Units', 'normalized',...
        'Position', [0.15 0.50 0.06 0.03],...
        'Min', 1,...
        'Max',  size(data.ct,2),...
        'Value', data.LateralOffset,...
        'SliderStep',[1/(size(data.ct,2)-1), 1/(size(data.ct,2)-1)],...
        'Callback', @ProfilesliderCallback,...
        'Enable',TmpStr,...
        'Visible',TmpStr);    
    
%% definition of callbacks
function ProfilesliderCallback(hObj,event)
      data=guidata(gcf);
      data.LateralOffset = round(get(hObj,'Value'));
      guidata(gcf,data);
      matRad_visCtDose;
end

 function displaypopupCallback(hObj,event)
      data=guidata(gcf);
      data.TypeOfPlot = get(hObj,'Value');
     
      guidata(gcf,data);
      matRad_visCtDose;
 
 end





    function dosepopupCallback(hObj,event)
      data=guidata(gcf);
      data.SelectedDisplayOption = data.fName{get(hObj,'Value'),1};
      guidata(gcf,data);
      matRad_visCtDose;
 
     end


     function doseColorwashCheckboxCallback(hObj,event)
        data = guidata(gcf);
        data.doseColorwashCheckboxValue = get(hObj,'Value');
        guidata(gcf,data);
        matRad_visCtDose;
     end
 
     function doseIsoCheckboxCallback(hObj,event)
        data = guidata(gcf);
        data.doseIsoCheckboxValue = get(hObj,'Value');
        guidata(gcf,data);
        matRad_visCtDose;
     end
 
    function ctCheckboxCallback(hObj,event)
        data = guidata(gcf);
        data.ctCheckboxValue = get(hObj,'Value');
        guidata(gcf,data);
        matRad_visCtDose;
    end

    function contourCheckboxCallback(hObj,event)
        data = guidata(gcf);
        data.contourCheckboxValue = get(hObj,'Value');
        guidata(gcf,data);
        matRad_visCtDose;
    end

    function planepopupCallback(hObj,event)
        data = guidata(gcf);
        plane = get(hObj,'Value');
        if data.plane ~= plane;
            data.plane = plane;
            data.slice = round(data.pln.isoCenter(data.plane)/data.pln.resolution(data.plane));
            if data.plane == 1
                axis([1 size(data.ct,3) 1 size(data.ct,1)]);
            elseif data.plane == 2        
                axis([1 size(data.ct,3) 1 size(data.ct,2)]);
            elseif data.plane == 3
                axis([1 size(data.ct,2) 1 size(data.ct,1)]);
            end
            guidata(gcf,data);
            matRad_visCtDose;
        end
    end

    function sliderCallback(hObj, event)
        data = guidata(gcf);
        data.slice = get(hObj, 'Value');
        data.slice = floor(data.slice);
        guidata(gcf,data);
        matRad_visCtDose;
    end



end
