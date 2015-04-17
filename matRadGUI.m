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

% Last Modified by GUIDE v2.5 17-Apr-2015 09:03:44

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
%% access variables from workspace
x=evalin('base','variableUsedinGUI');

if ~isempty(varargin)
    if ~isempty(varargin{1,1})
        handles.optResult = varargin{1,1};
        handles.State = 2;
    else
        handles.State = 1;
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
    

    % set slice slider
    handles.plane = get(handles.popupPlane,'value');
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

%% generate steering file
handles.stf = matRad_generateStf(handles.ct,handles.cst,handles.pln);
h=waitbar(0,'dose calculation ... ');
if strcmp(handles.pln.radiationMode,'photons')
    handles.dij = matRad_calcPhotonDose(handles.ct,handles.stf,handles.pln,handles.cst,0);
elseif strcmp(handles.pln.radiationMode,'protons') || strcmp(handles.pln.radiationMode,'carbon')
    handles.dij = matRad_calcParticleDose(handles.ct,handles.stf,handles.pln,handles.cst,0);
end
close(h);
handles.State = 2;
guidata(hObject,handles);
UpdatePlot(handles);


function UpdatePlot(handles)
cla;
plane=get(handles.popupPlane,'Value');
slice = get(handles.sliderSlice,'Value');
CutOffLevel= 0.05;

%% plot ct
if ~isempty(handles.ct.cube) && get(handles.popupTypeOfPlot,'Value')==1

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
if ~isempty(data.optResult) && data.TypeOfPlot ==1 && get(handles.radiobtnDose,'Value')

end
















%% plot VOIs
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

%% Set axis labels
if  plane == 1% Axial plane
    if ~isempty(handles.pln)
        set(gca,'XTick',0:50/handles.ct.resolution(1):1000)
        set(gca,'YTick',0:50/handles.ct.resolution(2):1000)
        set(gca,'XTickLabel',0:50:1000*handles.ct.resolution(1))
        set(gca,'YTickLabel',0:50:1000*handles.ct.resolution(2))
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
        xlabel('z [mm]','FontSize',16)
        ylabel('y [mm]','FontSize',16)
        title('Sagital plane','FontSize',15);
    else
        xlabel('z [voxels]','FontSize',16)
        ylabel('y [voxels]','FontSize',16)
        title('Sagital plane','FontSize',15);
    end
elseif plane == 3 % Coronal plane
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

















switch handles.State
    
    case 1 
        
    case 2
            
            
          
            handles.optResult.physicalDose=reshape(handles.dij.physicalDose*ones(handles.dij.totalNumOfBixels,1),handles.dij.dimensions);
            mVolume = handles.optResult.physicalDose;
            mVolume(handles.optResult.physicalDose<CutOffLevel*max(handles.optResult.physicalDose(:)))=0;
            
            dose_rgb = mVolume./max(mVolume(:));

            % Save RGB indices for dose in zsliceÂ´s voxels.
            if plane == 1  % Coronal plane
                
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(slice,:,:))),jet);
            elseif plane == 2 % Sagital plane
                
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,slice,:))),jet);
            elseif plane == 3 % Axial plane
                
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,slice))),jet);
            end
            
            doseImageHandle=imshow(dose_rgb, 'Parent', handles.axesFig);
   
            % Make dose transparent
            if ~isempty(handles.ct.cube)
                %set(doseImageHandle,'AlphaData',.6);
                if plane == 1  % Coronal plane
                    set(doseImageHandle,'AlphaData', .6*double(squeeze(handles.optResult.physicalDose(slice,:,:))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
                elseif plane == 2 % Sagital plane
                    set(doseImageHandle,'AlphaData', .6*double(squeeze(handles.optResult.physicalDose(:,slice,:))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
                elseif plane == 3 % Axial plane
                        set(doseImageHandle,'AlphaData', .6*double(squeeze(handles.optResult.physicalDose(:,:,slice))>CutOffLevel*max(handles.optResult.physicalDose(:))  )  ) ;
                end

            end
        
    case 3
        
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


% --- Executes on button press in radiobtnDose.
function radiobtnDose_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnDose


% --- Executes on button press in radiobtnIsoDoseLines.
function radiobtnIsoDoseLines_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnIsoDoseLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnIsoDoseLines


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
