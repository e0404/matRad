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

% Last Modified by GUIDE v2.5 04-May-2015 18:14:49

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
%show logo
axes(handles.axesLogo)
h=imshow('matrad_hat40.png');


%% 
% handles.State=0   no data available
% handles.State=1   data available ready for dose calculation
% handles.State=2   data available and dose is calculated;ready for
%                   optimization
% handles.Sate=3    plan is optimized


% if plan is changed go back to state 1
% if table is changed including VOI Type  go back to state 1
% if table is changed but not the VOI Types  go back to state 2
handles.TableChanged = false;
handles.State = 0;
%B = evalin('base','dij');
%assignin('base','dij2',B);
try
    if ~isempty(evalin('base','ct'))
        ct = evalin('base','ct');
        handles.State = 1;
    end
catch
    
end

try
    if ~isempty(evalin('base','doseVis'))
        handles.State = 2;

    end
catch 
end

try
    if ~isempty(evalin('base','optResult'))
        handles.State = 3;
    end
catch
end
%set cst
try
    if ~isempty(evalin('base','cst'))
        setCstTable(handles,evalin('base','cst'));
    end
catch
end
%set plan
try
    if ~isempty(evalin('base','pln'))
            pln=evalin('base','pln');
            setPln(handles);       
    end
catch
end

if handles.State ==2 || handles.State ==3
    set(handles.popupDisplayOption,'String','Dose');
    handles.SelectedDisplayOption ='Dose';
    handles.SelectedDisplayOptionIdx=1;
else
    handles.optResult = [];
    set(handles.popupDisplayOption,'String','no option available');
    handles.SelectedDisplayOption='';
    handles.SelectedDisplayOptionIdx=1;
end

handles.SelectedBeam=1;
handles.plane = get(handles.popupPlane,'Value');

if handles.State >0
    % set slice slider
    set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
        'Value',round(pln.isoCenter(handles.plane)/ct.resolution(handles.plane)),...
         'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
end

handles.profileOffset = 0;

% Update handles structure
guidata(hObject, handles);
UpdateState(handles)
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

% delete existing workspace
try
    clear global 'optResult'
catch 
end

try
    clear global 'doseVis'
catch 
end

try
    clear global 'dij'
catch 
end


[FileName, FilePath] = uigetfile;
load([FilePath FileName]);

setCstTable(handles,cst);

handles.TableChanged = false;
handles.State = 1;
set(handles.popupTypeOfPlot,'Value',1);
assignin('base','ct',ct);
assignin('base','cst',cst);

getPln(handles);
pln=evalin('base','pln');
% set slice slider
handles.plane = get(handles.popupPlane,'value');
set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
    'Value',round(pln.isoCenter(handles.plane)/ct.resolution(handles.plane)),...
     'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
  

guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);





function editSAD_Callback(hObject, eventdata, handles)
% hObject    handle to editSAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSAD as text
%        str2double(get(hObject,'String')) returns contents of editSAD as a double
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end
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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

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
if handles.State>0
    handles.State=1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in btnCalcDose.
function btnCalcDose_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalcDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% get cst from table
if ~getCstTable(handles);
    return
end
%% read plan from gui and save it to workspace
getPln(handles);
%% generate steering file
stf = matRad_generateStf(evalin('base','ct'),...
                                 evalin('base','cst'),...
                                 evalin('base','pln'));
assignin('base','stf',stf);
if strcmp(evalin('base','pln.radiationMode'),'photons')
    dij = matRad_calcPhotonDose(evalin('base','ct'),stf,evalin('base','pln'),evalin('base','cst'),0);
elseif strcmp(evalin('base','pln.radiationMode'),'protons') || strcmp(evalin('base','pln.radiationMode'),'carbon')
    dij = matRad_calcParticleDose(evalin('base','ct'),stf,evalin('base','pln'),evalin('base','cst'),0);
end

doseVis = matRad_mxCalcDose(dij,ones(dij.totalNumOfBixels,1),evalin('base','cst'));
assignin('base','dij',dij);
assignin('base','doseVis',doseVis);
handles.State = 2;
handles.SelectedDisplayOptionIdx=1;
handles.SelectedDisplayOption='Dose';
handles.SelectedBeam=1;
handles.TableChanged = false;
guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);


function UpdatePlot(handles)
%%
cla(handles.axesFig,'reset');

if handles.State ==0
    return
elseif handles.State > 0
     ct=evalin('base','ct');
     cst=evalin('base','cst');
     pln=evalin('base','pln');
end

if handles.State == 2
      Result = evalin('base','doseVis');
elseif handles.State == 3
      Result = evalin('base','optResult');
end

if exist('Result')
    if ~isempty(Result) && ~isempty(ct.cube)
        if isfield(Result,'w')
            Result = rmfield(Result,'w');
        end
        if isfield(Result,'physicalDose')
            Result.Dose = Result.physicalDose;
            Result=rmfield(Result,'physicalDose');
        end
        if isfield(Result,'RBE')
            Result.RBETruncated10Perc = Result.RBE;
            Result.RBETruncated10Perc(Result.Dose<0.1*...
                max(Result.Dose(:))) = 0;
        end

        fName =fieldnames(Result);
        for i=1:size(fName,1)
            %second dimension indicates if it should be plotted later on
            if strcmp(fName{i,1},'RBETruncated10Perc') || strcmp(fName{i,1},'w')
                fName{i,2}=0;
            else
                fName{i,2}=1;
            end
            % determine units
            if strcmp(fName{i,1},'Dose')
                fName{i,3} = '[Gy]';
            elseif strcmp(fName{i,1},'alpha')
                fName{i,3} = '[Gy^{-1}]';
            elseif strcmp(fName{i,1},'beta')
                fName{i,3} = '[Gy^{-2}]';
            elseif strcmp(fName{i,1},'RBEWeightedDose')
                fName{i,3} = '[Gy(RBE)]';
            else
                fName{i,3} = '[a.u.]';
            end
            % Reshape dose to cube in voxel dimensions
            CurrentCube = getfield(Result,fName{i,1});
            if ~isempty(CurrentCube) && isequal(size(CurrentCube),size(ct.cube))
                Result = setfield(Result,fName{i,1},reshape(CurrentCube,size(ct.cube)));
            %try to reshape using voxelDimensions from pln struct    
            elseif ~isempty(Result) && isequal(size(CurrentCube),size(ct.cube))
                Result = setfield(Result,fName{i,1},reshape(CurrentCube,pln.voxelDimensions));
            elseif ~isempty(Result) && ~strcmp(fName{i,1},'w')
                error('Cannot reshape dose');   
            end
        end

    set(handles.popupDisplayOption,'String',fName(:,1));
    set(handles.popupDisplayOption,'Value',find(strcmp(handles.SelectedDisplayOption,fName(:,1))));

    end

end




%% set and get required variables

plane=get(handles.popupPlane,'Value');
slice = round(get(handles.sliderSlice,'Value'));
CutOffLevel= 0.03;

%% plot ct
 if ~isempty(evalin('base','ct')) && get(handles.popupTypeOfPlot,'Value')==1
    cla(handles.axesFig);
    if plane == 1 % Coronal plane
        ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube(slice,:,:))/max(ct.cube(:))),bone);
    elseif plane == 2 % sagittal plane
        ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube(:,slice,:))/max(ct.cube(:))),bone);
    elseif plane == 3 % Axial plane
        ct_rgb = ind2rgb(uint8(63*squeeze(ct.cube(:,:,slice))/max(ct.cube(:))),bone);
    end
    
    axes(handles.axesFig)
    ctImageHandle = image(ct_rgb);
end

%% plot dose cube
if handles.State >1 &&  get(handles.popupTypeOfPlot,'Value')== 1 ...
        && get(handles.radiobtnDose,'Value')

        mVolume = getfield(Result,handles.SelectedDisplayOption);
        % make sure to exploit full color range
        mVolume(Result.Dose<CutOffLevel*max(Result.Dose(:)))=0;

    %     %% dose colorwash
        if ~isempty(mVolume) && get(handles.radiobtnDose,'Value') && ~isvector(mVolume)

            dose_rgb = mVolume./max(mVolume(:));

            % Save RGB indices for dose in zsliceÂ´s voxels.
            if plane == 1  % Coronal plane
                 mSlice = squeeze(mVolume(slice,:,:));
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(slice,:,:))),jet);
            elseif plane == 2 % sagittal plane
                 mSlice = squeeze(mVolume(:,slice,:));
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,slice,:))),jet);
            elseif plane == 3 % Axial plane
                 mSlice = squeeze(mVolume(:,:,slice));
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,slice))),jet);
            end
            % Show dose      
            doseImageHandle = image('CData',dose_rgb,'Parent',handles.axesFig);
            % make dose transparent
            if ~isempty(ct.cube)
                if plane == 1  % Coronal plane
                    set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.Dose(slice,:,:))>CutOffLevel*max(Result.Dose(:))  )  ) ;
                elseif plane == 2 % sagittal plane
                    set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.Dose(:,slice,:))>CutOffLevel*max(Result.Dose(:))  )  ) ;
                elseif plane == 3 % Axial plane
                    if strcmp(get(handles.popupDisplayOption,'String'),'RBETruncated10Perc')
                        set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.Dose(:,:,slice))>0.1*max(Result.Dose(:))  )  ) ;
                    else
                        set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.Dose(:,:,slice))>CutOffLevel*max(Result.Dose(:))  )  ) ;
                    end
                end

            end

            % plot colorbar
            v=version;
            if str2num(v(1:3))>=8.5
                cBarHandel = colorbar(handles.axesFig,'colormap',jet,'FontSize',8,'yAxisLocation','right');
            else
                cBarHandel = colorbar('peer',handles.axesFig,'FontSize',8,'yAxisLocation','right');
            end
            Idx = find(strcmp(handles.SelectedDisplayOption,fName(:,1)));
            set(get(cBarHandel,'ylabel'),'String', [fName{Idx,1} ' ' fName{Idx,3} ],'fontsize',8);

            if isempty(strfind(handles.SelectedDisplayOption,'RBE'))
                set(cBarHandel,'YLim',[0 max(eps,max(mVolume(:)))]);
                caxis(handles.axesFig,[0, max(eps,max(mVolume(:)))])
            else
                set(cBarHandel,'YLim',[0 max(eps,max(mSlice(:)))]);
                caxis(handles.axesFig,[0, max(eps,max(mSlice(:)))])
            end



        end
        
    axes(handles.axesFig),hold on 

        %% dose iso dose lines
        if get(handles.radiobtnIsoDoseLines,'Value')
                colormap(jet)
                vLevels = round(linspace(0,ceil(max(mVolume(:))),6).*100)/100;
                if plane == 1  % Coronal plane
                    Slice=squeeze(mVolume(slice,:,:));
                    if sum(Slice(:))>1
                        [~,myContour] = contour(Slice);
                    end
                elseif plane == 2 % Sagittal plane
                    Slice=squeeze(mVolume(:,slice,:));
                    if sum(Slice(:))>1
                        [~,myContour] = contour(Slice);
                    end
                elseif plane == 3 % Axial plane
                    Slice=squeeze(mVolume(:,:,slice));
                    if sum(Slice(:))>1
                        hold on
                     [~,myContour] = contour(Slice);
                    end
                end

                 if sum(Slice(:))>1
                    caxis(handles.axesFig,[0, max(mVolume(:))]);
                    % turn off legend for this data set
                    hAnnotation = get(myContour,'Annotation');
                    hLegendEntry = get(hAnnotation','LegendInformation');
                    set(hLegendEntry,'IconDisplayStyle','off');
                    set(myContour,'LabelSpacing',100);
                    if get(handles.radiobtnIsoDoseLinesLabels,'Value') == 0
                        set(myContour,'ShowText','off')
                    else
                        set(myContour,'ShowText','on')
                    end
                 end
        end

end


%% plot VOIs

if get(handles.radiobtnContour,'Value') && get(handles.popupTypeOfPlot,'Value')==1 && handles.State>0
    colors = jet;
    hold on,
    colors = colors(round(linspace(1,63,size(cst,1))),:);
    mask = zeros(size(ct.cube)); % create zero cube with same dimeonsions like dose cube
    for s = 1:size(cst,1)
        if ~strcmp(cst{s,3},'IGNORED') %&& ~strcmp(data.cst{s,2},'DoseFalloff')
            mask(:) = 0;
            mask(cst{s,4}) = 1;
            if plane == 1 && sum(sum(mask(slice,:,:))) > 0
                contour(handles.axesFig,squeeze(mask(slice,:,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',cst{s,2});
            elseif plane == 2 && sum(sum(mask(:,slice,:))) > 0
                contour(handles.axesFig,squeeze(mask(:,slice,:)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',cst{s,2});
            elseif plane == 3 && sum(sum(mask(:,:,slice))) > 0
                contour(handles.axesFig,squeeze(mask(:,:,slice)),.5*[1 1],'Color',colors(s,:),'LineWidth',2,'DisplayName',cst{s,2});
            end
        end
    end
    warning('off','MATLAB:legend:PlotEmpty')
    myLegend = legend('show','location','NorthEast');
    set(myLegend,'FontSize',8);
    set(myLegend,'color','none');
    set(myLegend,'TextColor', [1 1 1]);
    legend boxoff
    warning('on','MATLAB:legend:PlotEmpty')
end

%% Set axis labels
if  plane == 3% Axial plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution(1):1000);
        set(handles.axesFig,'YTick',0:50/ct.resolution(2):1000);
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution(1));
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution(2));   
        xlabel('x [mm]','FontSize',8)
        ylabel('y [mm]','FontSize',8)
        title('axial plane','FontSize',8)
    else
        xlabel('x [voxels]','FontSize',8)
        ylabel('y [voxels]','FontSize',8)
        title('axial plane','FontSize',8)
    end
elseif plane == 2 % Sagittal plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution(3):1000)
        set(handles.axesFig,'YTick',0:50/ct.resolution(2):1000)
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution(3))
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution(2))
        xlabel('z [mm]','FontSize',8);
        ylabel('y [mm]','FontSize',8);
        title('sagittal plane','FontSize',8);
    else
        xlabel('z [voxels]','FontSize',8)
        ylabel('y [voxels]','FontSize',8)
        title('sagittal plane','FontSize',8);
    end
elseif plane == 1 % Coronal plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution(3):1000)
        set(handles.axesFig,'YTick',0:50/ct.resolution(1):1000)
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution(3))
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution(1))
        xlabel('z [mm]','FontSize',8)
        ylabel('x [mm]','FontSize',8)
        title('coronal plane','FontSize',8)
    else
        xlabel('z [voxels]','FontSize',8)
        ylabel('x [voxels]','FontSize',8)
        title('coronal plane','FontSize',8)
    end
end



%% profile plot
if get(handles.popupTypeOfPlot,'Value')==2 && exist('Result')
     
    % clear view and initialize some values
    cla(handles.axesFig,'reset')
    set(gca,'YDir','normal');
    ylabel('{\color{black}dose [Gy]}')
    cColor={'black','green','magenta','cyan','yellow','red','blue'};
    
    % Rotation around Z axis (table movement)
    rotMx_XY = [ cosd(pln.gantryAngles(handles.SelectedBeam)) sind(pln.gantryAngles(handles.SelectedBeam)) 0;
                 -sind(pln.gantryAngles(handles.SelectedBeam)) cosd(pln.gantryAngles(handles.SelectedBeam)) 0;
                      0                                         0                                              1];
    % Rotation around Y axis (Couch movement)
     rotMx_XZ = [cosd(pln.couchAngles(handles.SelectedBeam)) 0 sind(pln.couchAngles(handles.SelectedBeam));
                 0                                             1 0;
                 -sind(pln.couchAngles(handles.SelectedBeam)) 0 cosd(pln.couchAngles(handles.SelectedBeam))];
    
    if strcmp(handles.ProfileType,'longitudinal')
        sourcePointBEV = [handles.profileOffset -pln.SAD   0];
        targetPointBEV = [handles.profileOffset pln.SAD   0];
        sMargin = -1;
    elseif strcmp(handles.ProfileType,'lateral')
        sourcePointBEV = [-pln.SAD handles.profileOffset   0];
        targetPointBEV = [pln.SAD handles.profileOffset   0];
        sMargin = 30;
    end
    rotSourcePointBEV = sourcePointBEV*rotMx_XY*rotMx_XZ;
    rotTargetPointBEV = targetPointBEV*rotMx_XY*rotMx_XZ;
    [~,~,~,~,vis] = matRad_siddonRayTracer(pln.isoCenter,ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{ct.cube},true);
    ix = vis.ix;
    mPhysDose=getfield(Result,'Dose'); 
    vPhysDose = mPhysDose(ix);
    % plot physical dose
    vX=linspace(1,ct.resolution(1)*numel(vPhysDose),numel(vPhysDose));
    PlotHandles{1} = plot(handles.axesFig,vX,vPhysDose,'color',cColor{1,1},'LineWidth',3);grid on, hold on; 
    PlotHandles{1,2}='Dose';
    set(gca,'FontSize',8);
    % assess x - limits
    xLim  = find(vPhysDose);
    if ~isempty(xLim)
        xmin= xLim(1)*ct.resolution(1)-sMargin;
        xmax= xLim(end)*ct.resolution(1)+sMargin;
    else
        vLim = axis;
        xmin = vLim(1);
        xmax = vLim(2);
    end
    
    % plot counter
    Cnt=2;
    
    if isfield(Result,'RBE')
        
        %disbale specific plots
        %data.fName{6,2}=0;
        %data.fName{5,2}=0;
        %data.fName{2,2}=0;
        
        % generate two lines for ylabel
        StringYLabel1 = '\fontsize{8}{\color{red}RBE x dose [Gy(RBE)] \color{black}dose [Gy] ';
        StringYLabel2 = '';
        for i=1:1:size(fName,1)
            mCurrentCube = getfield(Result,fName{i,1});
            if ~isvector(mCurrentCube) && ~strcmp(fName{i,1},'RBEWeightedDose') ...
                    && ~strcmp(fName{i,1},'RBE') && ~strcmp(fName{i,1},'Dose')...
                    && fName{i,2}
                vProfile = mCurrentCube(ix);
                PlotHandles{Cnt,1} = plot(vX,vProfile,'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                PlotHandles{Cnt,2} =fName{i,1};
                if strcmp(fName{i,1},'effect')
                    unit = 'a.u.';
                elseif strcmp(fName{i,1},'alpha')
                    unit = 'Gy^{-1}';
                elseif strcmp(fName{i,1},'beta')
                    unit = 'Gy^{-2}';
                end
              
                StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' fName{i,1} ' [' unit ']'];
                
                Cnt = Cnt+1;
            end           
        end
        StringYLabel2 = [StringYLabel2 '}'];
        % plot always RBEWeightedDose against RBE
        mRBEWeightedDose=getfield(Result,'RBEWeightedDose');
        vBED =mRBEWeightedDose(ix);
  
        mRBE=getfield(Result,'RBE');
        vRBE =mRBE(ix);
        
        % plot biological dose against RBE
        [ax, PlotHandles{Cnt,1}, PlotHandles{Cnt+1,1}]=plotyy(vX,vBED,vX,vRBE,'plot');hold on;
        PlotHandles{Cnt,2}='RBEWeightedDose';
        PlotHandles{Cnt+1,2}='RBE';
         
        % set plotyy properties
        set(get(ax(2),'Ylabel'),'String','RBE [a.u.]','FontSize',8);       
        ylabel({StringYLabel1;StringYLabel2})
        set(PlotHandles{Cnt,1},'Linewidth',4,'color','r');
        set(PlotHandles{Cnt+1,1},'Linewidth',3,'color','b');
        set(ax(1),'ycolor','r')
        set(ax(2),'ycolor','b')
        set(ax,'FontSize',8);
        Cnt=Cnt+2;
       
    end
       
    
    % asses target coordinates 
    tmpPrior = intmax;
    tmpSize = 0;
    for i=1:size(cst,1)
        if strcmp(cst{i,3},'TARGET') && tmpPrior>=cst{i,5}.Priority && tmpSize<numel(cst{i,4})
           mTarget = unique(cst{i,4});
           tmpPrior=cst{i,5}.Priority;
           tmpSize=numel(cst{i,4});
           VOI = cst{i,2};
        end
    end
    
    str = sprintf('profile plot of zentral axis of %d beam gantry angle %d° couch angle %d°',...
        handles.SelectedBeam ,pln.gantryAngles(handles.SelectedBeam),pln.couchAngles(handles.SelectedBeam));
    title(str,'FontSize',8),grid on
    
    
    % plot target boundaries
    mTargetStack = zeros(size(ct.cube));
    mTargetStack(mTarget)=1;
    vProfile =mTargetStack(ix);
    vRay = find(vProfile)*ct.resolution(2);
    
    PlotHandles{Cnt,2} =[VOI ' boundary'];
    vLim = axis;
    if ~isempty(vRay)
        PlotHandles{Cnt,1}=plot([vRay(1) vRay(1)],[0 vLim(4)],'--','Linewidth',3,'color','k');hold on
        plot([vRay(end) vRay(end)], [0 vLim(4)],'--','Linewidth',3,'color','k');hold on
    else
        PlotHandles{Cnt,1} =[];
    end
    
   Lines = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),1);
   Labels = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),2);
   h=legend([Lines{:}],Labels{:});
    set(h,'FontSize',8);
    % set axis limits
    if pln.bioOptimization == 0 || ~isfield(Result,'RBE')
        xlim([xmin xmax]);
     
    else
        xlim(ax(1),[xmin xmax]);
        xlim(ax(2),[xmin xmax]);
    end
    xlabel('depth [cm]','FontSize',8);
   
end

%% display command window output to textbox
jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
jCmdWin = jDesktop.getClient('Command Window');
jTextArea = jCmdWin.getComponent(0).getViewport.getView;
cwText = char(jTextArea.getText);
set(handles.listBoxCmd,'String',cwText);
numLines = size(get(handles.listBoxCmd,'String'),1);
if numLines>7
    set(handles.listBoxCmd, 'ListBoxTop', numLines-7);
end














% --- Executes on selection change in popupPlane.
function popupPlane_Callback(hObject, eventdata, handles)
% hObject    handle to popupPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupPlane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlane

contents = cellstr(get(hObject,'String'));
val = contents{get(hObject,'Value')};

% set slice slider
handles.plane = get(handles.popupPlane,'value');
% set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
%     'Value',round(pln.isoCenter(handles.plane)/ct.resolution(handles.plane)),...
%      'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
%  
try
    ct = evalin('base', 'ct');
    pln = evalin('base', 'pln');

    set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
       'Value',round(pln.isoCenter(handles.plane)/ct.resolution(handles.plane)),...
       'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
catch
end
        
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
set(handles.figure1,'Pointer','watch');
Param.numOfIter = str2num(get(handles.editNumIter,'String'));
Param.prec = str2num(get(handles.txtPrecisionOutput,'String'));
optResult = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),evalin('base','pln'),Param);
assignin('base','optResult',optResult);
set(handles.figure1,'Pointer','arrow');
handles.State=3;
handles.SelectedDisplayOptionIdx=1;
handles.SelectedDisplayOption='Dose';
handles.SelectedBeam=1;
guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);


% --- Executes on selection change in popupTypeOfPlot.
function popupTypeOfPlot_Callback(hObject, eventdata, handles)
% hObject    handle to popupTypeOfPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupTypeOfPlot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupTypeOfPlot
if get(hObject,'Value') ==1   % intensity plot
    set(handles.sliderBeamSelection,'Enable','off')
    set(handles.sliderOffset,'Enable','off')
    set(handles.popupDisplayOption,'Enable','on')
    set(handles.btnProfileType,'Enable','off');
    set(handles.popupPlane,'Enable','on');
    set(handles.radiobtnContour,'Enable','on');
    set(handles.radiobtnDose,'Enable','on');
    set(handles.radiobtnIsoDoseLines,'Enable','on');
    set(handles.radiobtnIsoDoseLinesLabels,'Enable','on');
    set(handles.sliderSlice,'Enable','on');
    
elseif get(hObject,'Value') ==2 % profile plot
    
    if handles.State >0
        if length(str2double(strsplit(get(handles.editGantryAngle,'String'),' ')))>1
            set(handles.sliderBeamSelection,'Enable','on');
            handles.SelectedBeam = 1;

                pln = evalin('base','pln');
                set(handles.sliderBeamSelection,'Min',handles.SelectedBeam,'Max',pln.numOfBeams,...
                    'Value',handles.SelectedBeam,...
                    'SliderStep',[1/(pln.numOfBeams-1) 1/(pln.numOfBeams-1)],...
                    'Enable','on');

        else
            handles.SelectedBeam=1;
        end
    
        handles.profileOffset=get(handles.sliderOffset,'Value');

        vMinMax = [-100 100];
        vRange = sum(abs(vMinMax));

        ct = evalin('base','ct');
        if strcmp(get(handles.btnProfileType,'String'),'lateral')
            SliderStep = vRange/ct.resolution(1);       
        else
            SliderStep = vRange/ct.resolution(2);  
        end
        
        set(handles.sliderOffset,'Min',vMinMax(1),'Max',vMinMax(2),...
                    'Value',handles.profileOffset,...
                    'SliderStep',[1/SliderStep 1/SliderStep],...
                    'Enable','on');
    end
    
    
    set(handles.popupDisplayOption,'Enable','off');
    set(handles.btnProfileType,'Enable','on');
    set(handles.popupPlane,'Enable','off');
    set(handles.radiobtnContour,'Enable','off');
    set(handles.radiobtnDose,'Enable','off');
    set(handles.radiobtnIsoDoseLines,'Enable','off');
    set(handles.sliderSlice,'Enable','off');
    set(handles.radiobtnIsoDoseLinesLabels,'Enable','off');
    
    
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


function setCstTable(handles,cst)


columnname = {'VOI','VOI Type','Overlap','Obj Func','Penalty','Dose','Parameters'};
columnformat = {cst(:,2)',{'OAR','TARGET'},'numeric',...
       {'square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
       'numeric','numeric','numeric'};
   
numOfObjectives = 0;
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        numOfObjectives = numOfObjectives + numel(cst{i,6});
    end
end

%dimArr = [size(cst,1)-sum(strcmp([cst(:,3)],'IGNORED')) size(columnname,2)];
dimArr = [numOfObjectives size(columnname,2)];
data = cell(dimArr);
Counter = 1;
for i = 1:size(cst,1)
   
   if strcmp(cst(i,3),'IGNORED')~=1
       for j=1:size(cst{i,6},1)
       %VOI
       data{Counter,1}=cst{i,2};
       %VOI Type
       data{Counter,2}=cst{i,3};
       %Priority
       data{Counter,3}=cst{i,5}.Priority;
       %Objective Function
       objFunc = cst{i,6}(j).type;
       data{Counter,4}=objFunc;
       %penalty
       data{Counter,5}=cst{i,6}(j).parameter(1);
       switch objFunc
           case 'EUD'
                data{Counter,7}=cst{i,6}(j).parameter(1);
                data{Counter,6}='';
           case 'mean'
               data{Counter,6}='';
               data{Counter,7}='';
           case {'square underdosing','square overdosing','square deviation'}
               %Dose
               data{Counter,6}=cst{i,6}(j).parameter(2);
               data{Counter,7}='';
       end
   
       Counter = Counter +1;
       end
   end
   
end
set(handles.uiTable,'ColumnName',columnname);
set(handles.uiTable,'ColumnFormat',columnformat);
set(handles.uiTable,'ColumnEditable',[true true true true true true true]);
% set(handles.uiTable,'ColumnWidth',{'auto','auto','auto','auto','auto','auto','auto'});
set(handles.uiTable,'Data',data);

function Flag=getCstTable (handles)

data = get(handles.uiTable,'Data');
OldCst = evalin('base','cst');
NewCst=[];
Cnt=1;
FlagValidParameters = true;
%% generate new cst from GUI
for i = 1:size(OldCst,1)
    CntObjF = 1;
    FlagFound = false;
    for j = 1:size(data,1)
        
        if strcmp(OldCst{i,2},data{j,1})
            FlagFound = true;
            
            if CntObjF == 1
                %VOI
                if isempty(data{j,1}) || ~isempty(findstr(data{j,1}, 'Select'))
                    FlagValidParameters=false;
                else
                    NewCst{Cnt,1}=data{j,1};
                end
                %VOI Type
                if isempty(data{j,2})|| ~isempty(findstr(data{j,2}, 'Select'))
                    FlagValidParameters=false;
                else
                    NewCst{Cnt,2}=data{j,2};
                end
                %Priority
                if isempty(data{j,3})
                    FlagValidParameters=false;
                else
                    NewCst{Cnt,3}=data{j,3};
                end
            end
            
            %Obj Func
            if isempty(data{j,4}) ||~isempty(findstr(data{j,4}, 'Select'))
               FlagValidParameters=false;
            else
                 NewCst{Cnt,4}(CntObjF,1).type = data{j,4};
            end
             %Penalty
            if isempty(data{j,5})
               FlagValidParameters=false;
            else
                 NewCst{Cnt,4}(CntObjF,1).parameter(1,1) = data{j,5};
            end
             
            %get exponent
            if FlagValidParameters
            
                if strcmp(NewCst{Cnt,4}(CntObjF,1).type,'EUD')
                    if isempty(data{j,7})
                       FlagValidParameters=false;
                    else
                        NewCst{Cnt,4}(CntObjF,1).exponent = data{j,7};
                    end
                end

                %get dose
                if sum(strcmp(NewCst{Cnt,4}(CntObjF,1).type,{'EUD','mean'})) == 0
                    % read dose
                    if isempty(data{j,6})
                       FlagValidParameters=false;
                    else
                          if length(NewCst{Cnt,4}(CntObjF,1).parameter)==1
                              NewCst{Cnt,4}(CntObjF,1).parameter(1,2)=1;
                          end
                          if iscellstr(data{j,6})
                              NewCst{Cnt,4}(CntObjF,1).parameter(1,2) = str2num(data{j,6});
                          else
                              NewCst{Cnt,4}(CntObjF,1).parameter(1,2) = (data{j,6});
                          end
                    end
                end
            end
            CntObjF = CntObjF+1; 
            
        end

    end
    
    if FlagFound == true
       Cnt = Cnt +1;
    end
            
end

if FlagValidParameters
    
       for m=1:size(OldCst,1)
           VOIexist   = OldCst(m,2);
           boolChanged = false;

           for n = 1:size(NewCst,1)

               VOIGUI = NewCst(n,1);

               if strcmp(VOIexist,VOIGUI)
                  % overite existing objectives
                   boolChanged = true;
                   OldCst(m,6)=NewCst(n,4);
                   OldCst(m,3)=NewCst(n,2);
                   OldCst{m,5}.Priority = NewCst{n,3};
               end 
           end

           if ~boolChanged
               OldCst{m,6}=[];
               OldCst{m,5}.Priority=nan;
           end

       end
       assignin('base','cst',OldCst);
       Flag = true;
       
else       
  warndlg('not all values are set - cannot start dose calculation'); 
  Flag = false;
end
%% replace old cst by new cst values and check for deleted objectives

% --- Executes on button press in btnuiTableAdd.
function btnuiTableAdd_Callback(hObject, eventdata, handles)
% hObject    handle to btnuiTableAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.uiTable, 'data');
sEnd = size(data,1);
data{sEnd+1,1} = 'Select VOI';
data{sEnd+1,2} = 'Select VOI Type';
data{sEnd+1,3} = 2;
data{sEnd+1,4} = 'Select obj func';
set(handles.uiTable,'data',data);

handles.State=1;
guidata(hObject,handles);
UpdateState(handles);


% --- Executes on button press in btnuiTableDel.
function btnuiTableDel_Callback(hObject, eventdata, handles)
% hObject    handle to btnuiTableDel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.uiTable, 'data');
rows = get(handles.uiTable,'UserData');
mask = (1:size(data,1))';
mask(rows)=[];
data=data(mask,:);
set(handles.uiTable,'data',data);
handles.State=1;
guidata(hObject,handles);
UpdateState(handles);

% --- Executes when selected cell(s) is changed in uiTable.
function uiTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uiTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
index = eventdata.Indices;
    if any(index)             %loop necessary to surpress unimportant errors.
        rows = index(:,1);
        set(hObject,'UserData',rows);
    end
    

% --- Executes when entered data in editable cell(s) in uiTable.
function uiTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uiTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% apply changes to the other cells
data = get(hObject,'Data');
if eventdata.Indices(2) == 1 || eventdata.Indices(2) == 2 ...
        || eventdata.Indices(2) == 3
    handles.State=1;
    handles.TableChanged = true;
else
    if handles.State ==3 && handles.TableChanged == false
        handles.State=2;
    end
end
%% check if input is a valid
if eventdata.Indices(2) == 3  || eventdata.Indices(2) == 5 || eventdata.Indices(2) == 6 || eventdata.Indices(2) == 7
    if CheckValidity(eventdata.NewData) ==false
            data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end
end

%% if objective function is set to mean --> set dose cell to empty
if eventdata.Indices(2)==4 && strcmp(eventdata.NewData,'mean') || ...
    eventdata.Indices(2)==4 && strcmp(eventdata.NewData,'EUD')
     data{eventdata.Indices(1),6}='';
end
%% if objective function is set to mean --> set EUD cell to empty
if eventdata.Indices(2)==4 && strcmp(eventdata.NewData,'mean')
    data{eventdata.Indices(1),7}='';
end


if eventdata.Indices(2) == 1 && eventdata.Indices(1) == size(data,1)
    for i = 1:size(data,1)
        if strcmp(eventdata.NewData,data{i,1})
           data{eventdata.Indices(1),2}=data{i,2};
           data{eventdata.Indices(1),3}=data{i,3};
        end
    end
    
end

%%
for i=1:size(data,1)
    if i~=eventdata.Indices(1) && strcmp(data(i,1),data(eventdata.Indices(1)))
        %set VOI type and priority
        data{i,2} = data{eventdata.Indices(1),2};
        data{i,3} = data{eventdata.Indices(1),3};
    end
end


%% check if editing current cell makes sense
%EUD column was edited
if eventdata.Indices(2) == 7
    if CheckValidity(eventdata.NewData) ==false
            data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end 
    % check if obj func is set to EUD otherwise reject this change
    if ~strcmp(data{eventdata.Indices(1),4},'EUD') 
        data{eventdata.Indices(1),eventdata.Indices(2)} = '';
    end
end

%Dose column was edited
if eventdata.Indices(2) == 6
    if CheckValidity(eventdata.NewData) ==false
            data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end 
    % check if obj func is set to EUD otherwise reject this change
    if sum(strcmp(data{eventdata.Indices(1),4},{'EUD','mean'}))>0
        data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end
end

if isnan(eventdata.NewData)
    data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
end

set(handles.txtInfo,'String','plan changed');
set(handles.uiTable,'data',data)
guidata(hObject, handles);
UpdateState(handles);



function FlagValidity = CheckValidity(Val) 
      
FlagValidity = true;

if  isempty(Val)
       warndlg('Input not a number !');
       FlagValidity = false;        
end

if Val < 0 
  warndlg('Input not a positive number !');
  FlagValidity = false;  
end

% if FlagValidity ==false
%       set(hObj,'BackgroundColor','r');
% else
%       set(hObj,'BackgroundColor',[0.94 0.94 0.94]);
% end
      
function UpdateState(handles)
 
 switch handles.State
     
     case 0
      
      set(handles.txtInfo,'String','no data loaded');
      set(handles.btnCalcDose,'Enable','off');
      set(handles.btnOptimize ,'Enable','off');
      set(handles.btnDVH,'Enable','off');
      
     case 1
      set(handles.txtInfo,'String','ready for dose calculation');
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','off');  
      set(handles.btnDVH,'Enable','off');
         
     case 2
      set(handles.txtInfo,'String','ready for optimization');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.btnDVH,'Enable','off');
     
     case 3
      set(handles.txtInfo,'String','plan is optimized');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.btnDVH,'Enable','on');
 end

 
 
 
function setPln(handles)
pln=evalin('base','pln');
set(handles.editBixelWidth,'String',num2str(pln.bixelWidth));
set(handles.editSAD,'String',num2str(pln.SAD));
set(handles.editFraction,'String',num2str(pln.numOfFractions));
set(handles.editGantryAngle,'String',num2str((pln.gantryAngles)));
set(handles.editCouchAngle,'String',num2str((pln.couchAngles)));
set(handles.popupRadMode,'Value',find(strcmp(get(handles.popupRadMode,'String'),pln.radiationMode)));
set(handles.radbtnBioOpt,'Value',pln.bioOptimization);
 
     
function getPln(handles)


pln.SAD             = parseStringAsNum(get(handles.editSAD,'String')); %[mm]
pln.bixelWidth      = parseStringAsNum(get(handles.editBixelWidth,'String')); % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = parseStringAsNum(get(handles.editGantryAngle,'String')); % [°]
pln.couchAngles     = parseStringAsNum(get(handles.editCouchAngle,'String')); % [°]

if length(pln.gantryAngles) ~= length(pln.couchAngles) 
  warndlg('number of gantryAngles != number of couchAngles'); 
end

pln.numOfBeams      = numel(pln.gantryAngles);
ct=evalin('base','ct');
pln.numOfVoxels     = numel(ct.cube);
pln.voxelDimensions = size(ct.cube);
contents                    = get(handles.popupRadMode,'String'); 
pln.radiationMode   =  contents{get(handles.popupRadMode,'Value')}; % either photons / protons / carbon
pln.bioOptimization = logical(get(handles.radbtnBioOpt,'Value'));   % false indicates physical optimization and true indicates biological optimization
pln.numOfFractions  = parseStringAsNum(get(handles.editFraction,'String'));
pln.voxelDimensions = size(ct.cube);
pln.isoCenter       = matRad_getIsoCenter(evalin('base','cst'),ct,0);
assignin('base','pln',pln);

function Number = parseStringAsNum(string)
try
    CellGantryAngles    = strsplit(string,'');
    Number = str2num(CellGantryAngles{1,1});
catch
      warndlg('could not parse all plan parameters'); 
end

% --- Executes on button press in btnTableSave.
function btnTableSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnTableSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getCstTable(handles);
getPln(handles);



% --- Executes on button press in btnDVH.
function btnDVH_Callback(hObject, eventdata, handles)
% hObject    handle to btnDVH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
matRad_calcDVH(evalin('base','optResult'),evalin('base','cst'))


% --- Executes on selection change in listBoxCmd.
function listBoxCmd_Callback(hObject, eventdata, handles)
% hObject    handle to listBoxCmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listBoxCmd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listBoxCmd
numLines = size(get(hObject,'String'),1);
set(hObject, 'ListboxTop', numLines);

% --- Executes during object creation, after setting all properties.
function listBoxCmd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listBoxCmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobtnIsoDoseLinesLabels.
function radiobtnIsoDoseLinesLabels_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnIsoDoseLinesLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnIsoDoseLinesLabels
UpdatePlot(handles);


% --- Executes on slider movement.
function sliderOffset_Callback(hObject, eventdata, handles)
% hObject    handle to sliderOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.profileOffset = get(hObject,'Value');
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function sliderOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function About_Callback(hObject, eventdata, handles)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myicon = imread('mR.png');
h=msgbox({'http://e0404.github.io/matRad/' 'm.banger@dkfz-heidelberg.de'},'About','custom',myicon);




% --- Executes on button press in btnRefresh.
function btnRefresh_Callback(hObject, eventdata, handles)
% hObject    handle to btnRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setPln(handles)
setCstTable(handles,evalin('base','cst'));
UpdatePlot(handles);


% --------------------------------------------------------------------
function toolbarSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbarSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
btnTableSave_Callback(hObject, eventdata, handles);



function editNumIter_Callback(hObject, eventdata, handles)
% hObject    handle to editNumIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumIter as text
%        str2double(get(hObject,'String')) returns contents of editNumIter as a double
val=round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function editNumIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slicerPrecision_Callback(hObject, eventdata, handles)
% hObject    handle to slicerPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
set(hObject,'Value',round(val))
set(handles.txtPrecisionOutput,'String',['1e-' num2str(round(val))]);

% --- Executes during object creation, after setting all properties.
function slicerPrecision_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slicerPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
