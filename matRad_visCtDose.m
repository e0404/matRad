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
    data.CutOffLevel = 0.03;
    data.optResult = optResult;
    data.cst  = cst;
    data.pln  = pln;
    data.ct   = ct;
    if ~isempty(data.optResult)
        data.optResult = rmfield(data.optResult,'w');
        if isfield(data.optResult,'RBE')
            data.optResult.RBETruncated10Perc = data.optResult.RBE;
            data.optResult.RBETruncated10Perc(data.optResult.physicalDose<0.1*max(data.optResult.physicalDose(:))) = 0;
        end
        
        data.fName =fieldnames(data.optResult);
        for i=1:size(data.fName,1)
            %indicates if it should be plotted later on
            if strcmp(data.fName{i,1},'RBETruncated10Perc')
                data.fName{i,2}=0;
            else
                data.fName{i,2}=1;
            end
            % determine units
            if strcmp(data.fName{i,1},'physicalDose')
                data.fName{i,3} = '[Gy]';
            elseif strcmp(data.fName{i,1},'alpha')
                data.fName{i,3} = '[Gy^{-1}]';
            elseif strcmp(data.fName{i,1},'beta')
                data.fName{i,3} = '[Gy^{-2}]';
            elseif strcmp(data.fName{i,1},'RBEWeightedDose')
                data.fName{i,3} = '[Gy(RBE)]';
            else
                data.fName{i,3} = '[a.u.]';
            end
            % Reshape dose to cube in voxel dimensions
            CurrentCube = getfield(data.optResult,data.fName{i,1});
            if ~isempty(CurrentCube) && ~isempty(data.ct.cube) && isequal(size(CurrentCube),size(data.ct.cube))
                data.optResult = setfield(data.optResult,data.fName{i,1},reshape(CurrentCube,size(ct.cube)));
            %try to reshape using voxelDimensions from pln struct    
            elseif ~isempty(data.optResult) && ~isempty(data.pln) && isequal(size(CurrentCube),size(data.ct.cube))
                data.optResult = setfield(data.optResult,data.fName{i,1},reshape(CurrentCube,data.pln.voxelDimensions));
            elseif ~isempty(data.optResult) && ~strcmp(data.fName{i,1},'w')
                error('Cannot reshape dose');   
            end
        end
    end
    
    data.SelectedBeam = 1;
    data.TypeOfPlot = 1;
    data.SelectedDisplayOption = 'physicalDose';
    data.ProfileType='lateral';
    if ~isempty(data.optResult)
        data.doseColorwashCheckboxValue = 1;
        data.doseIsoCheckboxValue = 1;
    else
        data.doseColorwashCheckboxValue = 0;
        data.doseIsoCheckboxValue = 0;        
    end
    
    if ~isempty(data.ct.cube)
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
        data.slice = round(pln.isoCenter(data.plane)/ct.resolution(data.plane));
    else
        data.slice = slice;
    end
    
    if data.plane == 1
        data.axis = [1 size(data.ct.cube,1) 1 size(data.ct.cube,3)];
    elseif data.plane == 2        
        data.axis = [1 size(data.ct.cube,3) 1 size(data.ct.cube,2)];
    elseif data.plane == 3
        data.axis = [1 size(data.ct.cube,1) 1 size(data.ct.cube,2)];
    end
    
    % Open figure
    myWindow = figure('Name','matRad CT/dose/VOI bowser','NumberTitle','off','units','normalized','outerposition',[0 0 1 1],'ToolBar','figure');
    myAxes   = axes('Position', [0.35 0.1 0.55 0.8],'YDir','reverse');
    guidata(myWindow,data);
    
else
    
    myWindow = gcf;
    data = guidata(myWindow);
    
    vAxis = [1 size(data.ct.cube,1) 1 size(data.ct.cube,2)];
    data.axis = vAxis;
    
    clf;
    
    myAxes   = axes('Position', [0.35 0.1 0.55 0.8],'YDir','reverse');
        
end

    hold on


    %% ct
if ~isempty(data.ct.cube) && data.ctCheckboxValue && data.TypeOfPlot ==1
    
    if data.plane == 1 % Coronal plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct.cube(data.slice,:,:))/max(data.ct.cube(:))),bone);
    elseif data.plane == 2 % Sagital plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct.cube(:,data.slice,:))/max(data.ct.cube(:))),bone);
    elseif data.plane == 3 % Axial plane
        ct_rgb = ind2rgb(uint8(63*squeeze(data.ct.cube(:,:,data.slice))/max(data.ct.cube(:))),bone);
    end
    axes(myAxes)
    ctImageHandle = image(ct_rgb);

end

if ~isempty(data.optResult) && data.TypeOfPlot ==1
    
    mVolume = getfield(data.optResult,data.SelectedDisplayOption);
    % make sure to exploit full color range
    mVolume(data.optResult.physicalDose<data.CutOffLevel*max(data.optResult.physicalDose(:)))=0;
    
%     %% dose colorwash
    if ~isempty(mVolume) && data.doseColorwashCheckboxValue && ~isvector(mVolume)

        dose_rgb = mVolume./max(mVolume(:));

        % Save RGB indices for dose in zsliceÂ´s voxels.
        if data.plane == 1  % Coronal plane
            mSlice = squeeze(mVolume(data.slice,:,:));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(data.slice,:,:))),jet);
        elseif data.plane == 2 % Sagital plane
             mSlice = squeeze(mVolume(:,data.slice,:));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,data.slice,:))),jet);
        elseif data.plane == 3 % Axial plane
             mSlice =squeeze(mVolume(:,:,data.slice));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,data.slice))),jet);
        end
        % Show dose
        axes(myAxes);
        doseImageHandle = image(dose_rgb);
        
        % Make dose transparent
        if ~isempty(data.ct.cube)
            %set(doseImageHandle,'AlphaData',.6);
        if data.plane == 1  % Coronal plane
            set(doseImageHandle,'AlphaData',  .6*double(squeeze(data.optResult.physicalDose(data.slice,:,:))>data.CutOffLevel*max(data.optResult.physicalDose(:))  )  ) ;
        elseif data.plane == 2 % Sagital plane
            set(doseImageHandle,'AlphaData',  .6*double(squeeze(data.optResult.physicalDose(:,data.slice,:))>data.CutOffLevel*max(data.optResult.physicalDose(:))  )  ) ;
        elseif data.plane == 3 % Axial plane
            if strcmp(data.SelectedDisplayOption,'RBETruncated10Perc')
                set(doseImageHandle,'AlphaData',  .6*double(squeeze(data.optResult.physicalDose(:,:,data.slice))>0.1*max(data.optResult.physicalDose(:))  )  ) ;
            else
                set(doseImageHandle,'AlphaData',  .6*double(squeeze(data.optResult.physicalDose(:,:,data.slice))>data.CutOffLevel*max(data.optResult.physicalDose(:))  )  ) ;
            end
        end
        
        end

    end

    %% dose iso dose lines
    if ~isempty(mVolume) && data.TypeOfPlot ==1 && ~isvector(mVolume)
           vLevels = round(linspace(0,ceil(max(mVolume(:))),4).*100)/100;
        if data.plane == 1  % Coronal plane
            [~,myContour] = contour(myAxes,squeeze(mVolume(data.slice,:,:)),vLevels);
        elseif data.plane == 2 % Sagittal plane
            [~,myContour] = contour(myAxes,squeeze(mVolume(:,data.slice,:)),vLevels);
        elseif data.plane == 3 % Axial plane
             [~,myContour] = contour(myAxes,squeeze(mVolume(:,:,data.slice)),vLevels);
        end

        % turn off legend for this data set
        hAnnotation = get(myContour,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
        set(myContour,'LabelSpacing',100,'ShowText','on')
        if data.doseIsoCheckboxValue == 0
            set(myContour,'Visible','off')
        end
    end

    cBarHandel = colorbar('peer',myAxes);
    Idx = find(strcmp(data.SelectedDisplayOption,data.fName(:,1)));
    set(get(cBarHandel,'ylabel'),'String', [data.fName{Idx,1} ' in ' data.fName{Idx,3} ],'fontsize',16);
    set(cBarHandel,'yAxisLocation','right');
    set(cBarHandel,'FontSize',14);
    
    if isempty(strfind(data.SelectedDisplayOption,'RBE'))
        set(cBarHandel,'YLim',[0 max(mVolume(:))]);
        caxis([0, max(mVolume(:))])
    else
        set(cBarHandel,'YLim',[0 max(mSlice(:))]);
        caxis([0, max(mSlice(:))])
    end

    
end
    %% VOIs
if ~isempty(data.cst) && data.contourCheckboxValue && data.TypeOfPlot ==1

    colors = jet;
    colors = colors(round(linspace(1,63,size(data.cst,1))),:);

    mask = zeros(size(data.ct.cube)); % create zero cube with same dimeonsions like dose cube
    for s = 1:size(data.cst,1)
        if ~strcmp(data.cst{s,3},'IGNORED') %&& ~strcmp(data.cst{s,2},'DoseFalloff')
            mask(:) = 0;
            mask(data.cst{s,4}) = 1;
            if data.plane == 1 && sum(sum(mask(data.slice,:,:))) > 0
                contour(myAxes,squeeze(mask(data.slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            elseif data.plane == 2 && sum(sum(mask(:,data.slice,:))) > 0
                contour(myAxes,squeeze(mask(:,data.slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            elseif data.plane == 3 && sum(sum(mask(:,:,data.slice))) > 0
                contour(myAxes,squeeze(mask(:,:,data.slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',data.cst{s,2});
            end
        end
    end

    myLegend = legend('show','location','NorthEast');
    %set(myLegend,'Position',[.85 .5 .1 .2]);
    set(myLegend,'FontSize',20);
    
end
    
%% Set axis labels
if   data.plane == 3% Axial plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.ct.resolution(1):1000)
        set(gca,'YTick',0:50/data.ct.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*data.ct.resolution(1))
        set(gca,'YTickLabel',0:50:1000*data.ct.resolution(2))
        xlabel('x [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        title('Axial plane','FontSize',16)
    else
        xlabel('x [voxels]','FontSize',16)
        ylabel('y [voxels]','FontSize',16)
        title('Axial plane','FontSize',16)
    end
    data.axis = [1 size(data.ct.cube,1) 1 size(data.ct.cube,2)];
elseif data.plane == 2 % Sagittal plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.ct.resolution(3):1000)
        set(gca,'YTick',0:50/data.ct.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*data.ct.resolution(3))
        set(gca,'YTickLabel',0:50:1000*data.ct.resolution(2))
        xlabel('z [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        title('Sagital plane','FontSize',15);
    else
        xlabel('z [voxels]','FontSize',16)
        ylabel('y [voxels]','FontSize',16)
        title('Sagital plane','FontSize',15);
    end
    data.axis = [1 size(data.ct.cube,3) 1 size(data.ct.cube,2)];
elseif data.plane == 1 % Coronal plane
    if ~isempty(data.pln)
        set(gca,'XTick',0:50/data.ct.resolution(3):1000)
        set(gca,'YTick',0:50/data.ct.resolution(1):1000)
        set(gca,'XTickLabel',0:50:1000*data.ct.resolution(3))
        set(gca,'YTickLabel',0:50:1000*data.ct.resolution(1))
        xlabel('z [mm]','FontSize',16)
        ylabel('x [mm]','FontSize',16)
        title('Coronal plane','FontSize',16)
    else
        xlabel('z [voxels]','FontSize',16)
        ylabel('x [voxels]','FontSize',16)
        title('Coronal plane','FontSize',16)
    end
    data.axis = [1 size(data.ct.cube,3) 1 size(data.ct.cube,1)];
end

%axis equal;
axis(data.axis);
set(gca,'FontSize',16);    


if data.TypeOfPlot ==2 &&~isempty(data.optResult)
     
    % clear view and initialize some values
    clf;
    axes('Position', [0.35 0.1 0.55 0.8]);
    set(gca,'YDir','normal');
    ylabel('{\color{black}dose in Gy}')
    cColor={'black','green','magenta','cyan','yellow','red','blue'};
    sSmoothFactor = 2;
    
    % Rotation around Z axis (table movement)
    rotMx_XY = [ cosd(data.pln.gantryAngles(data.SelectedBeam)) sind(data.pln.gantryAngles(data.SelectedBeam)) 0;
                 -sind(data.pln.gantryAngles(data.SelectedBeam)) cosd(data.pln.gantryAngles(data.SelectedBeam)) 0;
                      0                                         0                                              1];
    % Rotation around Y axis (Couch movement)
     rotMx_XZ = [cosd(data.pln.couchAngles(data.SelectedBeam)) 0 sind(data.pln.couchAngles(data.SelectedBeam));
                 0                                             1 0;
                 -sind(data.pln.couchAngles(data.SelectedBeam)) 0 cosd(data.pln.couchAngles(data.SelectedBeam))];
    if strcmp(data.ProfileType,'lateral')
        sourcePointBEV = [0 -data.pln.SAD   0];
        targetPointBEV = [0 data.pln.SAD   0];
        sMargin = -1;
    elseif strcmp(data.ProfileType,'longitudinal')
        sourcePointBEV = [-data.pln.SAD 0   0];
        targetPointBEV = [data.pln.SAD 0  0];
        sMargin = 30;
    end
    rotSourcePointBEV = sourcePointBEV*rotMx_XY*rotMx_XZ;
    rotTargetPointBEV = targetPointBEV*rotMx_XY*rotMx_XZ;
    [~,~,~,~,vis] = matRad_siddonRayTracer(data.pln.isoCenter,data.ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{data.ct.cube},true);
    ix = vis.ix;
    mPhysDose=getfield(data.optResult,'physicalDose'); %#ok<*GFLD>
    vPhysDose = mPhysDose(ix);
    % plot physical dose
    vX=linspace(1,data.ct.resolution(1)*numel(vPhysDose),numel(vPhysDose));
    
    PlotHandles{1} = plot(vX,smooth(vPhysDose,sSmoothFactor),'color',cColor{1,1},'LineWidth',3); hold on; 
    PlotHandles{1,2}='physicalDose';
    set(gca,'FontSize',18);
    % assess x - limits
    xLim  = find(vPhysDose);
    if ~isempty(xLim)
        xmin= xLim(1)*data.ct.resolution(1)-sMargin;
        xmax= xLim(end)*data.ct.resolution(1)+sMargin;
    else
        vLim = axis;
        xmin = vLim(1);
        xmax = vLim(2);
    end
    
    % plot counter
    Cnt=2;
    
    if isfield(data.optResult,'RBE')
        
        %disbale specific plots
        %data.fName{6,2}=0;
        %data.fName{5,2}=0;
        %data.fName{2,2}=0;
        
        % generate two lines for ylabel
        StringYLabel1 = '\fontsize{18}{\color{red}RBE x dose [Gy(RBE)] \color{black}physicalDose [Gy] ';
        StringYLabel2 = '';
        for i=1:1:size(data.fName,1)
            mCurrentCube = getfield(data.optResult,data.fName{i,1});
            if ~isvector(mCurrentCube) && ~strcmp(data.fName{i,1},'RBEWeightedDose') ...
                    && ~strcmp(data.fName{i,1},'RBE') && ~strcmp(data.fName{i,1},'physicalDose')...
                    && data.fName{i,2}
                vProfile = mCurrentCube(ix);
                PlotHandles{Cnt,1} = plot(vX,smooth(vProfile,sSmoothFactor),'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                PlotHandles{Cnt,2} =data.fName{i,1};
                if strcmp(data.fName{i,1},'effect')
                    unit = 'a.u.';
                elseif strcmp(data.fName{i,1},'alpha')
                    unit = 'Gy^{-1}';
                elseif strcmp(data.fName{i,1},'beta')
                    unit = 'Gy^{-2}';
                end
              
                StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' data.fName{i,1} ' [' unit ']'];
                
                Cnt = Cnt+1;
            end           
        end
        StringYLabel2 = [StringYLabel2 '}'];
        % plot always RBEWeightedDose against RBE
        mRBEWeightedDose=getfield(data.optResult,'RBEWeightedDose');
        vBED =mRBEWeightedDose(ix);
  
        mRBE=getfield(data.optResult,'RBE');
        vRBE =mRBE(ix);
        
        % plot biological dose against RBE
        [ax, PlotHandles{Cnt,1}, PlotHandles{Cnt+1,1}]=plotyy(vX,smooth(vBED,sSmoothFactor),vX,smooth(vRBE,sSmoothFactor),'plot');hold on;
        PlotHandles{Cnt,2}='RBEWeightedDose';
        PlotHandles{Cnt+1,2}='RBE';
         
        % set plotyy properties
        set(get(ax(2),'Ylabel'),'String','RBE [a.u.]','FontSize',18);       
        ylabel({StringYLabel1;StringYLabel2})
        set(PlotHandles{Cnt,1},'Linewidth',4,'color','r');
        set(PlotHandles{Cnt+1,1},'Linewidth',3,'color','b');
        set(ax(1),'ycolor','r')
        set(ax(2),'ycolor','b')
        set(ax,'FontSize',18);
        Cnt=Cnt+2;
       
    end
       
    
    % asses target coordinates 
    tmpPrior = intmax;
    tmpSize = 0;
    for i=1:size(data.cst,1)
        if strcmp(data.cst{i,3},'TARGET') && tmpPrior>=data.cst{i,5}.Priority && tmpSize<numel(data.cst{i,4})
           mTarget = unique(data.cst{i,4});
           tmpPrior=data.cst{i,5}.Priority;
           tmpSize=numel(data.cst{i,4});
           VOI = data.cst{i,2};
        end
    end
    

    str = sprintf('profile plot of zentral axis of %d beam gantry angle %d° couch angle %d°',...
        data.SelectedBeam ,data.pln.gantryAngles(data.SelectedBeam),data.pln.couchAngles(data.SelectedBeam));
    title(str,'FontSize',16),grid on
    
    
    % plot target boundaries
    mTargetStack = zeros(size(data.ct.cube));
    mTargetStack(mTarget)=1;
    vProfile =mTargetStack(ix);
    vRay = find(vProfile)*data.ct.resolution(2);
    
    PlotHandles{Cnt,2} =[VOI ' boundary'];
    vLim = axis;
    if ~isempty(vRay)
        PlotHandles{Cnt,1}=plot([vRay(1) vRay(1)],[0 vLim(4)],'--','Linewidth',3,'color','k');hold on
        plot([vRay(end) vRay(end)], [0 vLim(4)],'--','Linewidth',3,'color','k');hold on
        xmax = vRay(end)+30;
    else
        PlotHandles{Cnt,1} =0;
    end
    legend([PlotHandles{:,1}],PlotHandles{:,2},'Location','NorthWest');
    
    % set axis limits
    if data.pln.bioOptimization == 0 || ~isfield(data.optResult,'RBE')
        xlim([xmin xmax]);
    else
        xlim(ax(1),[xmin xmax]);
        xlim(ax(2),[xmin xmax]);
    end
    xlabel('depth [cm]','FontSize',16);
   
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
        'Max',  size(data.ct.cube,data.plane),...
        'Value', data.slice,...
        'SliderStep',[1/(size(data.ct.cube,data.plane)-1) 1/(size(data.ct.cube,data.plane)-1)],...
        'Callback', @sliderCallback);
    
doseSetText = uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Display options:',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.03 0.43 0.109 0.03],...
        'FontSize',14);   
    
if isempty(data.optResult)
    strTmp ={'no options available'};
    idx =1;
else
    strTmp = data.fName(:,1);
    idx = find(strcmp(data.fName(:,1),data.SelectedDisplayOption));
end

if data.TypeOfPlot == 2
    isEnabled = 'off';
else
    isEnabled = 'on';
end

uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', strTmp ,...
        'FontSize', 10,...
        'Value',idx,...
        'Units', 'normalized',...
        'Position', [0.05 0.4 0.109 0.03],...
        'Callback', @dosepopupCallback,...
        'Enable',isEnabled);

uicontrol('Parent', gcf,...
        'Style', 'text',...
        'BackgroundColor', [0.8 0.8 0.8],...
        'String', 'Type of plot:',...
        'FontSize', 10,...
        'Units', 'normalized',...
        'Position', [0.025 0.53 0.109 0.03],...
        'FontSize',14);    
    
uicontrol('Parent', gcf,...
        'Style', 'popupmenu',...
        'String', {'Intensity Plot','Profile Plot'} ,...
        'FontSize', 10,...
        'Value',data.TypeOfPlot,...
        'Units', 'normalized',...
        'Position', [0.05 0.50 0.059 0.03],...
        'Callback', @displaypopupCallback);

if data.TypeOfPlot == 2
    
    uicontrol('Parent', gcf,...
            'Style', 'togglebutton',...
            'String', data.ProfileType,...
            'Units', 'normalized',...
            'Position', [0.12 0.50 0.04 0.03],...
            'Callback', @toggleProfile);  
    
    
    if data.pln.numOfBeams>1

        uicontrol('Parent', gcf,...
            'Style', 'text',...
            'BackgroundColor', [0.8 0.8 0.8],...
            'String', 'beam selection',...
            'FontSize', 10,...
            'Units', 'normalized',...
            'Position', [0.17 0.43 0.109 0.03],...
            'FontSize',14);   

        uicontrol('Parent', gcf,...
                'Style', 'slider',...
                'Units', 'normalized',...
                'Position', [0.2 0.40 0.059 0.03],...
                'Min', 1,...
                'Max',  data.pln.numOfBeams,...
                'Value', data.SelectedBeam,...
                'SliderStep',[1/(data.pln.numOfBeams-1) 1/(data.pln.numOfBeams-1)],...
                'Callback', @beamSliderCallback);     
        
    end 
    
    
end   
%% definition of callbacks
 function toggleProfile(hObj,event)
      data=guidata(gcf);
      
      if strcmp(data.ProfileType,'longitudinal')
        set(hObj,'String','lateral');
      elseif strcmp(data.ProfileType,'lateral')
        set(hObj,'String','longitudinal');
      end
      data.ProfileType = get(hObj,'String');
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
      if isfield(data,'fName')
        data.SelectedDisplayOption = data.fName{get(hObj,'Value'),1};
      end
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
            data.slice = round(data.pln.isoCenter(data.plane)/data.ct.resolution(data.plane));
            if data.plane == 1
                axis([1 size(data.ct.cube,3) 1 size(data.ct.cube,1)]);
            elseif data.plane == 2        
                axis([1 size(data.ct.cube,3) 1 size(data.ct.cube,2)]);
            elseif data.plane == 3
                axis([1 size(data.ct.cube,2) 1 size(data.ct.cube,1)]);
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

 function beamSliderCallback(hObj, event)
        data = guidata(gcf);
        data.SelectedBeam = round(get(hObj, 'Value'));
        guidata(gcf,data);
        matRad_visCtDose;
    end


 end


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 