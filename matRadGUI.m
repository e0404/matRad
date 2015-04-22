function varargout = matRadGUI(varargin)
% MATRADGUI MATLAB code for matRadGUI.fig
%      MATRADGUI, by itself, creates a new MATRADGUI or raises the existing
%      singleton*.
%
%      H = MATRADGUI returns the handle to a new MATRADGUI or the handle to
%      the existing singleton*.
%
%      MATRADGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRADGUI.M with the given input arguments.
%
%      MATRADGUI('Property','Value',...) creates a new MATRADGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matRadGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matRadGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matRadGUI

% Last Modified by GUIDE v2.5 22-Apr-2015 15:57:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matRadGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @matRadGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before matRadGUI is made visible.
function matRadGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRadGUI (see VARARGIN)

% Choose default command line output for matRadGUI
handles.output = hObject;
 
%% 
if ~isempty(varargin)
    if ~isempty(varargin{1,1})
        handles.optResult = varargin{1,1};
        handles.State = 2;
        set(handles.popupDisplayOption,'String','physicalDose');
        handles.SelectedDisplayOption ='physicalDose';
        handles.SelectedDisplayOptionIdx=1;
    else
        handles.optResult = [];
        handles.State = 1;
        set(handles.popupDisplayOption,'String','no option available');
        handles.SelectedDisplayOption='';
        handles.SelectedDisplayOptionIdx=1;
    end
    if ~isempty(varargin{1,2})
        handles.cst = varargin{1,2};
    end
    if ~isempty(varargin{1,3})
        handles.pln = varargin{1,3};
        
        set(handles.editBixelWidth,'String',num2str(handles.pln.bixelWidth));
        set(handles.editSAD,'String',num2str(handles.pln.SAD));
        set(handles.editFraction,'String',num2str(handles.pln.numOfFractions));
        set(handles.editGantryAngle,'String',num2str(handles.pln.gantryAngles));
        set(handles.editCouchAngle,'String',num2str(handles.pln.couchAngles));
        set(handles.popupRadMode,'Value',find(strcmp(get(handles.popupRadMode,'String'),handles.pln.radiationMode)));
        set(handles.radbtnBioOpt,'Value',handles.pln.bioOptimization);
        
    end
    if ~isempty(varargin{1,4})
        handles.ct = varargin{1,4};
    end
    
    handles.pln.isoCenter = matRad_getIsoCenter(handles.cst,handles.ct,0);
    handles.plane = get(handles.popupPlane,'value');
    
    %% set beam slider
    handles.SelectedBeam=1;
    
    
    set(handles.sliderBeamSelection,'Min',handles.SelectedBeam,'Max',handles.pln.numOfBeams,...
        'Value',handles.SelectedBeam,...
        'SliderStep',[1/handles.pln.numOfBeams-1 1/handles.pln.numOfBeams-1],...
        'Enable','off');

    % set slice slider
    set(handles.sliderSlice,'Min',1,'Max',size(handles.ct.cube,handles.plane),...
        'Value',round(handles.pln.isoCenter(handles.plane)/handles.ct.resolution(handles.plane)),...
         'SliderStep',[1/(size(handles.ct.cube,handles.plane)-1) 1/(size(handles.ct.cube,handles.plane)-1)]);

end
% Update handles structure
guidata(hObject, handles);
UpdatePlot(handles)
% UIWAIT makes matRadGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matRadGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName FilePath] = uigetfile;
load([FilePath FileName]);
handles.pln.isoCenter = matRad_getIsoCenter(cst,ct,0);
handles.ct =ct;
handles.cst =cst;
handles.State = 1;

% set slice slider
handles.plane = get(handles.popupPlane,'value');
set(handles.sliderSlice,'Min',1,'Max',size(handles.ct.cube,handles.plane),...
    'Value',round(handles.pln.isoCenter(handles.plane)/handles.ct.resolution(handles.plane)),...
     'SliderStep',[1/(size(handles.ct.cube,handles.plane)-1) 1/(size(handles.ct.cube,handles.plane)-1)]);
  

guidata(hObject,handles);
UpdatePlot(handles);



function editSAD_Callback(hObject, eventdata, handles)
% hObject    handle to editSAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSAD as text
%        str2double(get(hObject,'String')) returns contents of editSAD as a double


% --- Executes during object creation, after setting all properties.
function editSAD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function editBixelWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editBixelWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBixelWidth as text
%        str2double(get(hObject,'String')) returns contents of editBixelWidth as a double


% --- Executes during object creation, after setting all properties.
function editBixelWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBixelWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGantryAngle_Callback(hObject, eventdata, handles)
% hObject    handle to editGantryAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGantryAngle as text
%        str2double(get(hObject,'String')) returns contents of editGantryAngle as a double


% --- Executes during object creation, after setting all properties.
function editGantryAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGantryAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCouchAngle_Callback(hObject, eventdata, handles)
% hObject    handle to editCouchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCouchAngle as text
%        str2double(get(hObject,'String')) returns contents of editCouchAngle as a double


% --- Executes during object creation, after setting all properties.
function editCouchAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCouchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupRadMode.
function popupRadMode_Callback(hObject, eventdata, handles)
% hObject    handle to popupRadMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupRadMode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupRadMode


% --- Executes during object creation, after setting all properties.
function popupRadMode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupRadMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFraction_Callback(hObject, eventdata, handles)
% hObject    handle to editFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFraction as text
%        str2double(get(hObject,'String')) returns contents of editFraction as a double


% --- Executes during object creation, after setting all properties.
function editFraction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radbtnBioOpt.
function radbtnBioOpt_Callback(hObject, eventdata, handles)
% hObject    handle to radbtnBioOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radbtnBioOpt


% --- Executes on button press in btnCalcDose.
function btnCalcDose_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalcDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.pln.SAD             = str2num(get(handles.editSAD,'String')); %[mm]
handles.pln.bixelWidth      = str2num(get(handles.editBixelWidth,'String')); % [mm] / also corresponds to lateral spot spacing for particles
handles.pln.gantryAngles    = str2num(get(handles.editGantryAngle,'String')); % [°]
handles.pln.couchAngles     = str2num(get(handles.editCouchAngle,'String')); % [°]
handles.pln.numOfBeams      = numel(handles.pln.gantryAngles);
handles.pln.numOfVoxels     = numel(handles.ct.cube);
handles.pln.voxelDimensions = size(handles.ct.cube);
contents                    = get(handles.popupRadMode,'String'); 
handles.pln.radiationMode   =  contents{get(handles.popupRadMode,'Value')}% either photons / protons / carbon
handles.pln.bioOptimization = logical(get(handles.radbtnBioOpt,'Value'));   % false indicates physical optimization and true indicates biological optimization
handles.pln.numOfFractions  = str2num(get(handles.editFraction,'String'));
handles.pln.voxelDimensions = size(handles.ct.cube);
%% generate steering file
handles.stf = matRad_generateStf(handles.ct,handles.cst,handles.pln);
h=waitbar(0,'dose calculation ... ');
if strcmp(handles.pln.radiationMode,'photons')
    handles.optResult = matRad_calcPhotonDose(handles.ct,handles.stf,handles.pln,handles.cst,0);
elseif strcmp(handles.pln.radiationMode,'protons') || strcmp(handles.pln.radiationMode,'carbon')
    handles.optResult = matRad_calcParticleDose(handles.ct,handles.stf,handles.pln,handles.cst,0);
end
close(h);
handles.State = 2;
guidata(hObject,handles);
UpdatePlot(handles);


function UpdatePlot(handles)
%%
    if ~isempty(handles.optResult) && ~isempty(handles.ct.cube)
        handles.optResult = rmfield(handles.optResult,'w');
        if isfield(handles.optResult,'RBE')
            handles.optResult.RBETruncated10Perc = handles.optResult.RBE;
            handles.optResult.RBETruncated10Perc(handles.optResult.physicalDose<0.1*...
                max(handles.optResult.physicalDose(:))) = 0;
        end
        
        handles.fName =fieldnames(handles.optResult);
        for i=1:size(handles.fName,1)
            %second dimension indicates if it should be plotted later on
            if strcmp(handles.fName{i,1},'RBETruncated10Perc') || strcmp(handles.fName{i,1},'w')
                handles.fName{i,2}=0;
            else
                handles.fName{i,2}=1;
            end
            % determine units
            if strcmp(handles.fName{i,1},'physicalDose')
                handles.fName{i,3} = '[Gy]';
            elseif strcmp(handles.fName{i,1},'alpha')
                handles.fName{i,3} = '[Gy^{-1}]';
            elseif strcmp(handles.fName{i,1},'beta')
                handles.fName{i,3} = '[Gy^{-2}]';
            elseif strcmp(handles.fName{i,1},'RBEWeightedDose')
                handles.fName{i,3} = '[Gy(RBE)]';
            else
                handles.fName{i,3} = '[a.u.]';
            end
            % Reshape dose to cube in voxel dimensions
            CurrentCube = getfield(handles.optResult,handles.fName{i,1});
            if ~isempty(CurrentCube) && isequal(size(CurrentCube),size(handles.ct.cube))
                handles.optResult = setfield(handles.optResult,handles.fName{i,1},reshape(CurrentCube,size(handles.ct.cube)));
            %try to reshape using voxelDimensions from pln struct    
            elseif ~isempty(handles.optResult) && isequal(size(CurrentCube),size(handles.ct.cube))
                handles.optResult = setfield(handles.optResult,handles.fName{i,1},reshape(CurrentCube,handles.pln.voxelDimensions));
            elseif ~isempty(handles.optResult) && ~strcmp(handles.fName{i,1},'w')
                error('Cannot reshape dose');   
            end
        end
    end




set(handles.popupDisplayOption,'String',handles.fName(:,1));
set(handles.popupDisplayOption,'Value',handles.SelectedDisplayOptionIdx);


%% set and get required variables

plane=get(handles.popupPlane,'Value');
slice = round(get(handles.sliderSlice,'Value'));
CutOffLevel= 0.03;

%% plot ct
 if ~isempty(handles.ct.cube) && get(handles.popupTypeOfPlot,'Value')==1
    cla;
    if plane == 1 % Coronal plane
        ct_rgb = ind2rgb(uint8(63*squeeze(handles.ct.cube(slice,:,:))/max(handles.ct.cube(:))),bone);
        axis([1 size(handles.ct.cube,1) 1 size(handles.ct.cube,2)]);
    elseif plane == 2 % Sagital plane
        ct_rgb = ind2rgb(uint8(63*squeeze(handles.ct.cube(:,slice,:))/max(handles.ct.cube(:))),bone);
        axis([1 size(handles.ct.cube,3) 1 size(handles.ct.cube,2)]);
    elseif plane == 3 % Axial plane
        ct_rgb = ind2rgb(uint8(63*squeeze(handles.ct.cube(:,:,slice))/max(handles.ct.cube(:))),bone);
        axis([1 size(handles.ct.cube,1) 1 size(handles.ct.cube,3)]); 
    end

    imshow(ct_rgb, 'Parent', handles.axesFig),hold on;
    
end

%% plot dose cube
if ~isempty(handles.optResult) &&  get(handles.popupTypeOfPlot,'Value')== 1 ...
    && get(handles.radiobtnDose,'Value')

    mVolume = getfield(handles.optResult,handles.SelectedDisplayOption);
    % make sure to exploit full color range
    mVolume(handles.optResult.physicalDose<CutOffLevel*max(handles.optResult.physicalDose(:)))=0;
    
%     %% dose colorwash
    if ~isempty(mVolume) && get(handles.radiobtnDose,'Value') && ~isvector(mVolume)

        dose_rgb = mVolume./max(mVolume(:));
    
        % Save RGB indices for dose in zsliceÂ´s voxels.
        if plane == 1  % Coronal plane
             mSlice = squeeze(mVolume(slice,:,:));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(slice,:,:))),jet);
        elseif plane == 2 % Sagital plane
             mSlice = squeeze(mVolume(:,slice,:));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,slice,:))),jet);
        elseif plane == 3 % Axial plane
             mSlice = squeeze(mVolume(:,:,slice));
            dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,slice))),jet);
        end
        % Show dose
        axes(handles.axesFig);
        doseImageHandle = image(dose_rgb);
        % make dose transparent
        if ~isempty(handles.ct.cube)
          
            if plane == 1  % Coronal plane
                set(doseImageHandle,'AlphaData',  .6*double(squeeze(handles.optResult.physicalDose(slice,:,:))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
            elseif plane == 2 % Sagital plane
                set(doseImageHandle,'AlphaData',  .6*double(squeeze(handles.optResult.physicalDose(:,slice,:))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
            elseif plane == 3 % Axial plane
                if strcmp(get(handles.popupDisplayOption,'String'),'RBETruncated10Perc')
                    set(doseImageHandle,'AlphaData',  .6*double(squeeze(handles.optResult.physicalDose(:,:,slice))>0.1*max(handles.optResult.physicalDose(:))  )  ) ;
                else
                    set(doseImageHandle,'AlphaData',  .6*double(squeeze(handles.optResult.physicalDose(:,:,slice))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
                end
            end
        
        end

        % plot colorbar
        cBarHandel = colorbar('peer',handles.axesFig,'colormap',jet,'FontSize',14,'yAxisLocation','right');
        Idx = find(strcmp(handles.SelectedDisplayOption,handles.fName(:,1)));
        set(get(cBarHandel,'ylabel'),'String', [handles.fName{Idx,1} ' in ' handles.fName{Idx,3} ],'fontsize',16);

        if isempty(strfind(handles.SelectedDisplayOption,'RBE'))
            set(cBarHandel,'YLim',[0 max(mVolume(:))]);
            caxis([0, max(mVolume(:))])
        else
            set(cBarHandel,'YLim',[0 max(mSlice(:))]);
            caxis([0, max(mSlice(:))])
        end
        
        
        
    end
    
    
    %% dose iso dose lines
    if get(handles.radiobtnIsoDoseLines,'Value')
            colormap(jet)
           
            vLevels = round(linspace(0,ceil(max(mVolume(:))),6).*100)/100;
            if plane == 1  % Coronal plane
                Slice=squeeze(mVolume(slice,:,:));
                if sum(Slice(:))>1
                    [~,myContour] = contour(handles.axesFig,Slice,vLevels);
                end
            elseif plane == 2 % Sagittal plane
                Slice=squeeze(mVolume(:,slice,:));
                if sum(Slice(:))>1
                    [~,myContour] = contour(handles.axesFig,Slice,vLevels);
                end
            elseif plane == 3 % Axial plane
                Slice=squeeze(mVolume(:,:,slice));
                if sum(Slice(:))>1
                 [~,myContour] = contour(handles.axesFig,Slice,vLevels);
                end
            end
            
             if sum(Slice(:))>1
                caxis(handles.axesFig,[0, max(mVolume(:))]);
                % turn off legend for this data set
                hAnnotation = get(myContour,'Annotation');
                hLegendEntry = get(hAnnotation','LegendInformation');
                set(hLegendEntry,'IconDisplayStyle','off')
                set(myContour,'LabelSpacing',100,'ShowText','on')
             end
    end
    
     
end


%% plot VOIs

if get(handles.radiobtnContour,'Value') && get(handles.popupTypeOfPlot,'Value')==1
    colors = jet;
    colors = colors(round(linspace(1,63,size(handles.cst,1))),:);
    mask = zeros(size(handles.ct.cube)); % create zero cube with same dimeonsions like dose cube
    for s = 1:size(handles.cst,1)
        if ~strcmp(handles.cst{s,3},'IGNORED') %&& ~strcmp(data.cst{s,2},'DoseFalloff')
            mask(:) = 0;
            mask(handles.cst{s,4}) = 1;
            if plane == 1 && sum(sum(mask(slice,:,:))) > 0
                contour(handles.axesFig,squeeze(mask(slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',handles.cst{s,2});
            elseif plane == 2 && sum(sum(mask(:,slice,:))) > 0
                contour(handles.axesFig,squeeze(mask(:,slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',handles.cst{s,2});
            elseif plane == 3 && sum(sum(mask(:,:,slice))) > 0
                contour(handles.axesFig,squeeze(mask(:,:,slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',handles.cst{s,2});
            end
        end
    end

    myLegend = legend('show','location','NorthEast');
    set(myLegend,'FontSize',12);
end

%% Set axis labels
if  plane == 3% Axial plane
    if ~isempty(handles.pln)
        set(handles.axesFig,'XTick',0:50/handles.ct.resolution(1):1000);
        set(handles.axesFig,'YTick',0:50/handles.ct.resolution(2):1000);
        set(handles.axesFig,'XTickLabel',0:50:1000*handles.ct.resolution(1));
        set(gca,'YTickLabel',0:50:1000*handles.ct.resolution(2));   
        xlabel('x [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        title('Axial plane','FontSize',16)
    else
        xlabel('x [voxels]','FontSize',16)
        ylabel('y [voxels]','FontSize',16)
        title('Axial plane','FontSize',16)
    end
elseif plane == 2 % Sagittal plane
    if ~isempty(handles.pln)
        set(gca,'XTick',0:50/handles.ct.resolution(3):1000)
        set(gca,'YTick',0:50/handles.ct.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*handles.ct.resolution(3))
        set(gca,'YTickLabel',0:50:1000*handles.ct.resolution(2))
        xlabel('z [mm]','FontSize',16);
        ylabel('y [mm]','FontSize',16);
        title('Sagital plane','FontSize',15);
    else
        xlabel('z [voxels]','FontSize',16)
        ylabel('y [voxels]','FontSize',16)
        title('Sagital plane','FontSize',15);
    end
elseif plane == 1 % Coronal plane
    if ~isempty(handles.pln)
        set(gca,'XTick',0:50/handles.ct.resolution(3):1000)
        set(gca,'YTick',0:50/handles.ct.resolution(1):1000)
        set(gca,'XTickLabel',0:50:1000*handles.ct.resolution(3))
        set(gca,'YTickLabel',0:50:1000*handles.ct.resolution(1))
        xlabel('z [mm]','FontSize',16)
        ylabel('x [mm]','FontSize',16)
        title('Coronal plane','FontSize',16)
    else
        xlabel('z [voxels]','FontSize',16)
        ylabel('x [voxels]','FontSize',16)
        title('Coronal plane','FontSize',16)
    end
end

axis equal;
set(gca,'FontSize',16);  


%% profile plot
if get(handles.popupTypeOfPlot,'Value')==2 && ~isempty(handles.optResult)
     
    % clear view and initialize some values
    cla(handles.axesFig,'reset')
    set(gca,'YDir','normal');
    ylabel('{\color{black}dose in Gy}')
    cColor={'black','green','magenta','cyan','yellow','red','blue'};
    sSmoothFactor = 2;
    
    % Rotation around Z axis (table movement)
    rotMx_XY = [ cosd(handles.pln.gantryAngles(handles.SelectedBeam)) sind(handles.pln.gantryAngles(handles.SelectedBeam)) 0;
                 -sind(handles.pln.gantryAngles(handles.SelectedBeam)) cosd(handles.pln.gantryAngles(handles.SelectedBeam)) 0;
                      0                                         0                                              1];
    % Rotation around Y axis (Couch movement)
     rotMx_XZ = [cosd(handles.pln.couchAngles(handles.SelectedBeam)) 0 sind(handles.pln.couchAngles(handles.SelectedBeam));
                 0                                             1 0;
                 -sind(handles.pln.couchAngles(handles.SelectedBeam)) 0 cosd(handles.pln.couchAngles(handles.SelectedBeam))];
    
    if strcmp(handles.ProfileType,'longitudinal')
        sourcePointBEV = [0 -handles.pln.SAD   0];
        targetPointBEV = [0 handles.pln.SAD   0];
        sMargin = -1;
    elseif strcmp(handles.ProfileType,'lateral')
        sourcePointBEV = [-handles.pln.SAD 0   0];
        targetPointBEV = [handles.pln.SAD 0  0];
        sMargin = 30;
    end
    rotSourcePointBEV = sourcePointBEV*rotMx_XY*rotMx_XZ;
    rotTargetPointBEV = targetPointBEV*rotMx_XY*rotMx_XZ;
    [~,~,~,~,vis] = matRad_siddonRayTracer(handles.pln.isoCenter,handles.ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{handles.ct.cube},true);
    ix = vis.ix;
    mPhysDose=getfield(handles.optResult,'physicalDose'); 
    vPhysDose = mPhysDose(ix);
    % plot physical dose
    vX=linspace(1,handles.ct.resolution(1)*numel(vPhysDose),numel(vPhysDose));
    PlotHandles{1} = plot(handles.axesFig,vX,smooth(vPhysDose,sSmoothFactor),'color',cColor{1,1},'LineWidth',3);grid on, hold on; 
    PlotHandles{1,2}='physicalDose';
    set(gca,'FontSize',18);
    % assess x - limits
    xLim  = find(vPhysDose);
    if ~isempty(xLim)
        xmin= xLim(1)*handles.ct.resolution(1)-sMargin;
        xmax= xLim(end)*handles.ct.resolution(1)+sMargin;
    else
        vLim = axis;
        xmin = vLim(1);
        xmax = vLim(2);
    end
    
    % plot counter
    Cnt=2;
    
    if isfield(handles.optResult,'RBE')
        
        %disbale specific plots
        %data.fName{6,2}=0;
        %data.fName{5,2}=0;
        %data.fName{2,2}=0;
        
        % generate two lines for ylabel
        StringYLabel1 = '\fontsize{18}{\color{red}RBE x dose [Gy(RBE)] \color{black}physicalDose [Gy] ';
        StringYLabel2 = '';
        for i=1:1:size(handles.fName,1)
            mCurrentCube = getfield(handles.optResult,handles.fName{i,1});
            if ~isvector(mCurrentCube) && ~strcmp(handles.fName{i,1},'RBEWeightedDose') ...
                    && ~strcmp(handles.fName{i,1},'RBE') && ~strcmp(handles.fName{i,1},'physicalDose')...
                    && handles.fName{i,2}
                vProfile = mCurrentCube(ix);
                PlotHandles{Cnt,1} = plot(vX,smooth(vProfile,sSmoothFactor),'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                PlotHandles{Cnt,2} =handles.fName{i,1};
                if strcmp(handles.fName{i,1},'effect')
                    unit = 'a.u.';
                elseif strcmp(handles.fName{i,1},'alpha')
                    unit = 'Gy^{-1}';
                elseif strcmp(handles.fName{i,1},'beta')
                    unit = 'Gy^{-2}';
                end
              
                StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' handles.fName{i,1} ' [' unit ']'];
                
                Cnt = Cnt+1;
            end           
        end
        StringYLabel2 = [StringYLabel2 '}'];
        % plot always RBEWeightedDose against RBE
        mRBEWeightedDose=getfield(handles.optResult,'RBEWeightedDose');
        vBED =mRBEWeightedDose(ix);
  
        mRBE=getfield(handles.optResult,'RBE');
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
    for i=1:size(handles.cst,1)
        if strcmp(handles.cst{i,3},'TARGET') && tmpPrior>=handles.cst{i,5}.Priority && tmpSize<numel(handles.cst{i,4})
           mTarget = unique(handles.cst{i,4});
           tmpPrior=handles.cst{i,5}.Priority;
           tmpSize=numel(handles.cst{i,4});
           VOI = handles.cst{i,2};
        end
    end
    
    str = sprintf('profile plot of zentral axis of %d beam gantry angle %d° couch angle %d°',...
        handles.SelectedBeam ,handles.pln.gantryAngles(handles.SelectedBeam),handles.pln.couchAngles(handles.SelectedBeam));
    title(str,'FontSize',16),grid on
    
    
    % plot target boundaries
    mTargetStack = zeros(size(handles.ct.cube));
    mTargetStack(mTarget)=1;
    vProfile =mTargetStack(ix);
    vRay = find(vProfile)*handles.ct.resolution(2);
    
    PlotHandles{Cnt,2} =[VOI ' boundary'];
    vLim = axis;
    if ~isempty(vRay)
        PlotHandles{Cnt,1}=plot([vRay(1) vRay(1)],[0 vLim(4)],'--','Linewidth',3,'color','k');hold on
        plot([vRay(end) vRay(end)], [0 vLim(4)],'--','Linewidth',3,'color','k');hold on
    else
        PlotHandles{Cnt,1} =0;
    end
    legend([PlotHandles{:,1}],PlotHandles{:,2});
    
    % set axis limits
    if handles.pln.bioOptimization == 0 || ~isfield(handles.optResult,'RBE')
        xlim([xmin xmax]);
     
    else
        xlim(ax(1),[xmin xmax]);
        xlim(ax(2),[xmin xmax]);
    end
    xlabel('depth [cm]','FontSize',12);
   
end

















% --- Executes on selection change in popupPlane.
function popupPlane_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupPlane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlane
UpdatePlot(handles)

% --- Executes during object creation, after setting all properties.
function popupPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderSlice_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
UpdatePlot(handles)

% --- Executes during object creation, after setting all properties.
function sliderSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobtnContour.
function radiobtnContour_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnContour
UpdatePlot(handles)

% --- Executes on button press in radiobtnDose.
function radiobtnDose_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnDose
UpdatePlot(handles)

% --- Executes on button press in radiobtnIsoDoseLines.
function radiobtnIsoDoseLines_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnIsoDoseLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnIsoDoseLines
UpdatePlot(handles)

% --- Executes on button press in btnOptimize.
function btnOptimize_Callback(hObject, eventdata, handles)
% hObject    handle to btnOptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=waitbar(0,'optimizaiton ... ');
handles.optResult = matRad_inversePlanning(handles.dij,handles.cst,handles.pln);
handles.State = 3;
close(h)
guidata(hObject,handles);
UpdatePlot(handles);


% --- Executes on selection change in popupTypeOfPlot.
function popupTypeOfPlot_Callback(hObject, eventdata, handles)
% hObject    handle to popupTypeOfPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupTypeOfPlot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupTypeOfPlot
if get(hObject,'Value') ==1   % intensity plot
    set(handles.sliderBeamSelection,'Enable','off')
    set(handles.popupDisplayOption,'Enable','on')
    set(handles.btnProfileType,'Enable','off');
    set(handles.popupPlane,'Enable','on');
    set(handles.radiobtnContour,'Enable','on');
    set(handles.radiobtnDose,'Enable','on');
    set(handles.radiobtnIsoDoseLines,'Enable','on');
    set(handles.sliderSlice,'Enable','on');
    
elseif get(hObject,'Value') ==2 % profile plot
    
    if handles.pln.numOfBeams>1
        set(handles.sliderBeamSelection,'Enable','on')
    end
    set(handles.popupDisplayOption,'Enable','off')
    set(handles.btnProfileType,'Enable','on');
    set(handles.popupPlane,'Enable','off');
    set(handles.radiobtnContour,'Enable','off');
    set(handles.radiobtnDose,'Enable','off');
    set(handles.radiobtnIsoDoseLines,'Enable','off');
    set(handles.sliderSlice,'Enable','off');
    
    handles.SelectedBeam = get(handles.sliderBeamSelection,'Value');
    
    set(handles.btnProfileType,'Enable','on')
    if strcmp(get(handles.btnProfileType,'String'),'lateral')
        handles.ProfileType = 'longitudinal';
    else
        handles.ProfileType = 'lateral';
    end
    
end
guidata(hObject, handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function popupTypeOfPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupTypeOfPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupDisplayOption.
function popupDisplayOption_Callback(hObject, eventdata, handles)
% hObject    handle to popupDisplayOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupDisplayOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDisplayOption
content = get(hObject,'String');
handles.SelectedDisplayOption = content{get(hObject,'Value'),1};
handles.SelectedDisplayOptionIdx = get(hObject,'Value');
guidata(hObject, handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function popupDisplayOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupDisplayOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderBeamSelection_Callback(hObject, eventdata, handles)
% hObject    handle to sliderBeamSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SelectedBeam = get(hObject,'Value');
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function sliderBeamSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderBeamSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnProfileType.
function btnProfileType_Callback(hObject, eventdata, handles)
% hObject    handle to btnProfileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(get(hObject,'Enable') ,'on')
    if strcmp(handles.ProfileType,'lateral')
            handles.ProfileType = 'longitudinal';
             set(hObject,'String','lateral');
        else
            handles.ProfileType = 'lateral';
            set(hObject,'String','longitudinal');
    end
 guidata(hObject, handles);
 UpdatePlot(handles);
 
end
