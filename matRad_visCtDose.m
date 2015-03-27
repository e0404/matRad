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
        data.fName =fieldnames(data.optResult);
        for i=1:size(data.fName,1)
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
        data.SelectedDisplayOption = 2;
        data.TypeOfPlot = 1;
    else
        data.doseColorwashCheckboxValue = 0;
        data.doseIsoCheckboxValue = 0;
        data.SelectedDisplayOption = 2;
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
    
    data.profileY = NaN;
    
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
    mVolume = getfield(data.optResult,data.fName{data.SelectedDisplayOption});
%     %% dose colorwash
    if ~isempty(mVolume) && data.doseColorwashCheckboxValue

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
    if ~isempty(mVolume) && data.doseIsoCheckboxValue && data.TypeOfPlot ==1
         
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
    ccc={'k','g','r','b','m','y','c'};
    ymax=0;
    %it is possible to add a certain width to the profile plot
    delta =0; % make it bixel distance dependend

    
    mDose=getfield(data.optResult,'physicalDose');
    mActualSlice = mDose(:,:,data.slice);
    mRotActualSlice =imrotate(mActualSlice,data.pln.gantryAngles(1),'crop');
    vW =ones(size(mRotActualSlice,2),1);
    
    vProjected =vW'*mRotActualSlice;
    %find first and last nonzero element of projected data 
    %which can be seen as determinng the beam width
    idx = find(vProjected);
    idxCentAxis = round((idx(end)+idx(1))/2);
    
    if isnan(data.profileY)
        data.profileY = idxCentAxis;
    else
        idxCentAxis=data.profileY;
    end
    
    % plot the physical dose
    vY=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
    vY(isnan(vY))=0;
    vY_avg=mean(vY,2);
    vX=linspace(1,data.pln.resolution(1)*numel(vY_avg),numel(vY_avg));
    h1=plot(vX,vY_avg,'color',ccc{1,1},'LineWidth',3),hold on; 
    
    % assess x and y axis limits
    xLim  = find(vY_avg);
    xmin= xLim(1)*data.pln.resolution(1)-20;
    xmax= xLim(end)*data.pln.resolution(1)+20;
    if max(vY_avg(:))>ymax
             ymax=max(vY_avg(:));
    end
    
    
    if data.pln.bioOptimization == 1 && strcmp(data.pln.radiationMode,'carbon')
        
        mEffect=getfield(data.optResult,'effect');
        mActualSlice = mEffect(:,:,data.slice);
        mRotActualSlice =imrotate(mActualSlice,data.pln.gantryAngles(1),'crop');
        mEffect=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mEffect(isnan(mEffect))=0;
        vEffect=mean(mEffect,2);
        h2=plot(vX,vEffect,'color',ccc{1,2},'LineWidth',3),hold on; 
        
        mAlpha=getfield(data.optResult,'alpha');
        mActualSlice = mAlpha(:,:,data.slice);
        mRotActualSlice =imrotate(mActualSlice,data.pln.gantryAngles(1),'crop');
        mAlpha=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mAlpha(isnan(mAlpha))=0;
        vAlpha=mean(mAlpha,2);
        h3=plot(vX,vAlpha,'color',ccc{1,7},'LineWidth',3),hold on; 
      
        mBioDose=getfield(data.optResult,'RBEWeightedDose');
        mActualSlice = mBioDose(:,:,data.slice);
        mRotActualSlice =imrotate(mActualSlice,data.pln.gantryAngles(1),'crop');
        mBioDose=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mBioDose(isnan(mBioDose))=0;
        vBioDose=mean(mBioDose,2);
        
        mRBE=getfield(data.optResult,'RBE');
        mActualSlice = mRBE(:,:,data.slice);
        mRotActualSlice =imrotate(mActualSlice,data.pln.gantryAngles(1),'crop');
        mRBE=mRotActualSlice(:,idxCentAxis-delta:idxCentAxis+delta);
        mRBE(isnan(mRBE))=0;
        vRBE=mean(mRBE,2);
        
        % plot biological dose against RBE
        [ax,h4,h5]=plotyy(vX,vBioDose,vX,vRBE,'plot'),hold on;
      
        
        set(get(ax(2),'Ylabel'),'String','RBE','FontSize',16);
        set(get(ax(1),'Ylabel'),'String','RBE x dose','FontSize',16);
        set(h4,'Linewidth',4,'color',ccc{1,3});
        set(h5,'Linewidth',3,'color',ccc{1,4});
        set(ax(1),'ycolor','r')
        set(ax(2),'ycolor','b')
        set(ax,'FontSize',16);
      
        
         if max(vBioDose(:))>ymax
             ymax=max(vBioDose(:));
         end
    end
       
    
    % asses the prescripted dose and target coordinates 
    % todo: ptv, ctv, gtv are all labeld as target
    sPrescrpDose=0;
    for i=1:size(data.cst,1)
        if strcmp(data.cst{i,3},'TARGET')==1
           mTarget = unique(data.cst{i,8});
           sPrescrpDose = data.cst{i,4};
        end
        
    end
    
    % plot target boundaries
    mTargetStack = zeros(size(data.ct));
    mTargetStack(mTarget)=1;
    mTargetSlice = mTargetStack(:,:,data.slice);
    mRotTargetSlice =imrotate(mTargetSlice,data.pln.gantryAngles(1),'crop');
    vRay = find(mRotTargetSlice(:,idxCentAxis))*data.pln.resolution(2);
    if ~isempty(vRay)
        h6=plot([vRay(1) vRay(1)],[0 ymax],'--','Linewidth',2,'color','k'),hold on
        plot([vRay(end) vRay(end)], [0 ymax],'--','Linewidth',2,'color','k'),hold on
        xmax = vRay(end)+30;
    else
        h6 =0;
    end
    
    % plot prescription
    h7=plot([0 size(data.ct,1)*data.pln.resolution(1)],[sPrescrpDose sPrescrpDose],'--','Linewidth',2,'color','m')
    
    str = sprintf('profile plot of zentral axis of first beam at %d° in slice %d / %d ',data.pln.gantryAngles(1),data.profileY*data.pln.resolution(2),size(data.ct,2)*data.pln.resolution(2));
    title(str,'FontSize',14),grid on
    axis auto
    if data.pln.bioOptimization == 1 && strcmp(data.pln.radiationMode,'carbon')
        legend([h1;h2;h3;h4;h5;h6;h7],'Physical Dose','Effect','Alpha','Biological Dose','RBE','target boundary','prescription');      
    else
        legend([h1;h6;h7],'Physical Dose','target boundary','prescription');
    end
    
    % set axis limits
    if data.pln.bioOptimization == 0 || ~strcmp(data.pln.radiationMode,'carbon')
        xlim([xmin xmax]);
    else
        xlim(ax(1),[xmin xmax]);
        xlim(ax(2),[xmin xmax]);
    end
    if ymax<0.5
        ymax = 0.5;
    end
    % set axis labels and ticks
    ymax= 0:sPrescrpDose/8:sPrescrpDose+sPrescrpDose/8;
    yVal=round(ymax*100)/100;
    set(gca,'YTick',yVal);
    set(gca,'YTickLabel',yVal);
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
    strTmp = data.fName;
end

dosePopup = uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', strTmp ,...
        'FontSize', 10,...
        'Value',data.SelectedDisplayOption,...
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
        'Value', data.profileY,...
        'SliderStep',[1/(size(data.ct,2)-1), 1/(size(data.ct,2)-1)],...
        'Callback', @ProfilesliderCallback,...
        'Enable',TmpStr,...
        'Visible',TmpStr);    
    
%% definition of callbacks
function ProfilesliderCallback(hObj,event)
      data=guidata(gcf);
      data.profileY = round(get(hObj,'Value'));
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
      data.SelectedDisplayOption = get(hObj,'Value');
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
