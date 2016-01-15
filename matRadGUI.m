function varargout = matRadGUI(varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad GUI
%
% call
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

% Begin initialization code - DO NOT EDIT
% set platform specific look and feel
if ispc
    lf = 'com.sun.java.swing.plaf.windows.WindowsLookAndFeel';
elseif isunix
    lf = 'com.jgoodies.looks.plastic.Plastic3DLookAndFeel';
elseif ismac
    lf = 'com.apple.laf.AquaLookAndFeel';
end
javax.swing.UIManager.setLookAndFeel(lf);

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
function matRadGUI_OpeningFcn(hObject, ~, handles, varargin) 
%#ok<*DEFNU> 
%#ok<*AGROW>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRadGUI (see VARARGIN)

% Choose default command line output for matRadGUI
handles.output = hObject;
%show matrad logo
axes(handles.axesLogo)
[im, ~, alpha] = imread(['dicomImport' filesep 'matrad_logo.png']);
f = image(im);
axis equal off
set(f, 'AlphaData', alpha);
% show dkfz logo
axes(handles.axesDKFZ)
[im, ~, alpha] = imread(['dicomImport' filesep 'DKFZ_logo.png']);
f = image(im);
axis equal off;
set(f, 'AlphaData', alpha);


%initialize maximum dose for visualization to Zero
handles.maxDoseVal     = 0;
handles.IsoDose.RefVal = 0;
handles.IsoDose.Levels = 0;
%seach for availabes machines
handles.Modalities = {'photons','protons','carbon'};
for i = 1:length(handles.Modalities)
    pattern = [handles.Modalities{1,i} '_*'];
    Files = dir(pattern);
    if isdeployed
        Files = [Files, dir([ctfroot filesep 'matRad' filesep pattern])];
    end
    for j = 1:length(Files)
        if ~isempty(Files)
            MachineName = Files(j).name(numel(handles.Modalities{1,i})+2:end-4);
            if isfield(handles,'Machines')
                if sum(strcmp(handles.Machines,MachineName)) == 0
                  handles.Machines{size(handles.Machines,1)+1} = MachineName;
                end
            else
                handles.Machines = cell(1);
                handles.Machines{1} = MachineName;
            end
        end
    end
end
set(handles.popUpMachine,'String',handles.Machines);


vChar = get(handles.editGantryAngle,'String');
if strcmp(vChar(1,1),'0') && length(vChar)==6
    set(handles.editGantryAngle,'String','0');
end
vChar = get(handles.editCouchAngle,'String');
if strcmp(vChar(1,1),'0') && length(vChar)==3
    set(handles.editCouchAngle,'String','0')
end
%% 
% handles.State=0   no data available
% handles.State=1   ct cst and pln available; ready for dose calculation
% handles.State=2   ct cst and pln available and dij matric(s) are calculated;
%                   ready for optimization
% handles.State=3   plan is optimized


% if plan is changed go back to state 1
% if table VOI Type or Priorities changed go to state 1
% if objective parameters changed go back to state 2

handles.TableChanged = false;
handles.State = 0;

%% parse variables from base workspace
AllVarNames = evalin('base','who');

try
    if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
        ct  = evalin('base','ct');
        cst = evalin('base','cst');
        setCstTable(handles,cst);
        handles.State = 1;
    elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find cst file');
    elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find ct file');
    end
catch  
   handles = showError(handles,'GUI OpeningFunc: Could not load ct and cst file');
end

%set plan if available - if not create one
try 
     if ismember('pln',AllVarNames) && handles.State > 0 
          setPln(handles); 
     elseif handles.State > 0 
          getPln(handles);
     end
catch
       handles.State = 0;
       handles = showError(handles,'GUI OpeningFunc: Could not set or get pln');
end

% check for dij structure
if ismember('dij',AllVarNames)
    handles.State = 2;
end

% check for optimized results
if ismember('resultGUI',AllVarNames)
    handles.State = 3;
end

% set some default values
if handles.State == 2 || handles.State == 3
    set(handles.popupDisplayOption,'String','physicalDose');
    handles.SelectedDisplayOption ='physicalDose';
    handles.SelectedDisplayOptionIdx=1;
else
    handles.resultGUI = [];
    set(handles.popupDisplayOption,'String','no option available');
    handles.SelectedDisplayOption='';
    handles.SelectedDisplayOptionIdx=1;
end

%per default the first beam is selected for the profile plot
handles.SelectedBeam = 1;
handles.plane = get(handles.popupPlane,'Value');
handles.DijCalcWarning = false;

    % set slice slider
if handles.State > 0
        set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
            'Value',ceil(size(ct.cube,handles.plane)/2),...
            'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);      
end

% Update handles structure
handles.profileOffset = 0;
guidata(hObject, handles);
UpdateState(handles)
UpdatePlot(handles)




% --- Outputs from this function are returned to the command line.
function varargout = matRadGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% set focus on error dialog
if isfield(handles,'ErrorDlg')
    figure(handles.ErrorDlg)
end

% --- Executes on button press in btnLoadMat.
function btnLoadMat_Callback(hObject, ~, handles)
% hObject    handle to btnLoadMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% delete existing workspace - parse variables from base workspace
AllVarNames = evalin('base','who');
RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
for i = 1:length(RefVarNames)  
    if sum(ismember(AllVarNames,RefVarNames{i}))>0
        evalin('base',['clear ', RefVarNames{i}]);
    end
end

% read new data
handles.State = 0;
try 
    [FileName, FilePath] = uigetfile;
    load([FilePath FileName]);
catch
    handles = showWarning(handles,'LoadMatFileFnc: Could not load *.mat file');
    guidata(hObject,handles);
    UpdatePlot(handles);
    UpdateState(handles);
    return
end

try
    setCstTable(handles,cst);
    handles.TableChanged = false;
    set(handles.popupTypeOfPlot,'Value',1);

    assignin('base','ct',ct);
    assignin('base','cst',cst);
catch
    handles = showError(handles,'LoadMatFileFnc: Could not load selected data');
end

try
    if exist('pln','var')
        % assess plan from loaded *.mat file
        assignin('base','pln',pln);
        setPln(handles);
    else
        % assess plan variable from GUI
        getPln(handles);
        setPln(handles);
    end
    handles.State = 1;
catch
    handles.State = 0;
end


% check if a optimized plan was loaded
if exist('stf','var')  && exist('dij','var')
    assignin('base','stf',stf);
    assignin('base','dij',dij);
    handles.State = 2;
end

if exist('resultGUI','var')
    assignin('base','resultGUI',resultGUI);
    handles.State = 3;
    handles.SelectedDisplayOption ='physicalDose';
end

% set slice slider
handles.plane = get(handles.popupPlane,'value');
if handles.State >0
     set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
            'Value',round(size(ct.cube,handles.plane)/2),...
            'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
end

guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);


% --- Executes on button press in btnLoadDicom.
function btnLoadDicom_Callback(hObject, ~, handles)
% hObject    handle to btnLoadDicom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    % delete existing workspace - parse variables from base workspace
    AllVarNames = evalin('base','who');
    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
    for i = 1:length(RefVarNames)  
        if sum(ismember(AllVarNames,RefVarNames{i}))>0
            evalin('base',['clear ', RefVarNames{i}]);
        end
    end
    handles.State = 0;
    if ~isdeployed
        addpath([pwd filesep 'dicomImport']);
    end
    matRad_importDicomGUI;
 
catch
   handles = showError(handles,'DicomImport: Could not import data'); 
end
UpdateState(handles);
guidata(hObject,handles);

function editBixelWidth_Callback(hObject, ~, handles)
% hObject    handle to editBixelWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBixelWidth as text
%        str2double(get(hObject,'String')) returns contents of editBixelWidth as a double
getPln(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function editBixelWidth_CreateFcn(hObject, ~, ~)
% hObject    handle to editBixelWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editGantryAngle_Callback(hObject, ~, handles)
% hObject    handle to editGantryAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGantryAngle as text
%        str2double(get(hObject,'String')) returns contents of editGantryAngle as a double
getPln(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function editGantryAngle_CreateFcn(hObject, ~, ~)
% hObject    handle to editGantryAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCouchAngle_Callback(hObject, ~, handles)
% hObject    handle to editCouchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCouchAngle as text
%        str2double(get(hObject,'String')) returns contents of editCouchAngle as a double
getPln(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function editCouchAngle_CreateFcn(hObject, ~, ~)
% hObject    handle to editCouchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupRadMode.
function popupRadMode_Callback(hObject, ~, handles)
% hObject    handle to popupRadMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkRadiationComposition(handles);
contents = cellstr(get(hObject,'String')); 
RadIdentifier = contents{get(hObject,'Value')};

switch RadIdentifier
    case 'photons'
        set(handles.radbtnBioOpt,'Value',0);
        set(handles.radbtnBioOpt,'Enable','off');
        set(handles.btnTypBioOpt,'Enable','off');
        
        set(handles.btnRunSequencing,'Enable','on');
        set(handles.btnRunDAO,'Enable','on');
        set(handles.txtSequencing,'Enable','on');
        set(handles.editSequencingLevel,'Enable','on');
        
    case 'protons'
        set(handles.radbtnBioOpt,'Value',0);
        set(handles.radbtnBioOpt,'Enable','off');
        set(handles.btnTypBioOpt,'Enable','off');
        
        set(handles.btnRunSequencing,'Enable','off');
        set(handles.btnRunDAO,'Enable','off');
        set(handles.txtSequencing,'Enable','off');
        set(handles.editSequencingLevel,'Enable','off');
        
    case 'carbon'
        set(handles.radbtnBioOpt,'Value',1);
        set(handles.radbtnBioOpt,'Enable','on');
        set(handles.btnTypBioOpt,'Enable','on');
        
        set(handles.btnRunSequencing,'Enable','off');
        set(handles.btnRunDAO,'Enable','off');
        set(handles.txtSequencing,'Enable','off');
        set(handles.editSequencingLevel,'Enable','off');
end

if handles.State > 0
    pln = evalin('base','pln');
    if handles.State > 0 && ~strcmp(contents(get(hObject,'Value')),pln.radiationMode)
        handles.State = 1;
        UpdateState(handles);
        guidata(hObject,handles);
    end
   getPln(handles);
end



% --- Executes during object creation, after setting all properties.
function popupRadMode_CreateFcn(hObject, ~, ~)
% hObject    handle to popupRadMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFraction_Callback(hObject, ~, handles)
% hObject    handle to editFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFraction as text
%        str2double(get(hObject,'String')) returns contents of editFraction as a double
getPln(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function editFraction_CreateFcn(hObject, ~, ~) 
% hObject    handle to editFraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radbtnBioOpt.
function radbtnBioOpt_Callback(hObject, ~, handles)
% hObject    handle to radbtnBioOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radbtnBioOpt
getPln(handles);
if get(hObject,'Value')
    set(handles.btnTypBioOpt,'Enable','on');
else
    set(handles.btnTypBioOpt,'Enable','off');
end

if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in btnCalcDose.
function btnCalcDose_Callback(hObject, ~, handles)
% hObject    handle to btnCalcDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% http://stackoverflow.com/questions/24703962/trigger-celleditcallback-before-button-callback
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/332613
% wait some time until the CallEditCallback is finished 
% Callback triggers the cst saving mechanism from GUI
try
    % indicate that matRad is busy
    % change mouse pointer to hour glass 
    Figures = findobj('type','figure');
    set(Figures, 'pointer', 'watch'); 
    drawnow;
    % disable all active objects
    InterfaceObj = findobj(Figures,'Enable','on');
    set(InterfaceObj,'Enable','off');
    
    pause(0.1);
    uiTable_CellEditCallback(hObject,[],handles);
    pause(0.3);

    %% get cst from table
    if ~getCstTable(handles);
        return
    end
    % read plan from gui and save it to workspace
    % gets also IsoCenter from GUI if checkbox is not checked
    getPln(handles);

    % get default iso center as center of gravity of all targets if not
    % already defined
    pln = evalin('base','pln');

    if length(pln.gantryAngles) ~= length(pln.couchAngles) 
        handles = showWarning(handles,warndlg('number of gantryAngles != number of couchAngles')); 
    end
    %%
    if ~checkRadiationComposition(handles);
        fileName = [pln.radiationMode '_' pln.machine];
        handles = showError(handles,errordlg(['Could not find the following machine file: ' fileName ]));
        guidata(hObject,handles);
        return;
    end

    %% check if isocenter is already set
    if ~isfield(pln,'isoCenter')
        handles = showWarning(handles,warning('no iso center set - using center of gravity based on structures defined as TARGET'));
        pln.isoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
        assignin('base','pln',pln);
    elseif ~get(handles.checkIsoCenter,'Value') 
        pln.isoCenter = str2num(get(handles.editIsoCenter,'String'));
    end

catch
   handles = showError(handles,'CalcDoseCallback: Error in preprocessing step.'); 
   guidata(hObject,handles);
   return;
end

% generate steering file
try 
    stf = matRad_generateStf(evalin('base','ct'),...
                                     evalin('base','cst'),...
                                     evalin('base','pln'));
    assignin('base','stf',stf);
catch
   handles = showError(handles,'CalcDoseCallback: Error in steering file generation'); 
   guidata(hObject,handles);
   return;
end

% carry out dose calculation
try
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(evalin('base','ct'),stf,pln,evalin('base','cst'),0);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(evalin('base','ct'),stf,pln,evalin('base','cst'),0);
    end

    % assign results to base worksapce
    assignin('base','dij',dij);
    handles.State = 2;
    handles.TableChanged = false;
    UpdatePlot(handles);
    UpdateState(handles);
    guidata(hObject,handles);
catch
    handles = showError(handles,'CalcDoseCallback: Error in dose calculation'); 
    guidata(hObject,handles);
    return;
end

% change state from busy to normal
set(Figures, 'pointer', 'arrow');
set(InterfaceObj,'Enable','on');

guidata(hObject,handles);  



%% plots ct and distributions
function UpdatePlot(handles)

defaultFontSize = 8;
cla(handles.axesFig,'reset');

if handles.State == 0
    return
elseif handles.State > 0
     ct  = evalin('base','ct');
     cst = evalin('base','cst');
     pln = evalin('base','pln');
end

%% state 3 indicates that a optimization has been performed
if handles.State > 2
      Result = evalin('base','resultGUI');
end

if exist('Result','var')
    if ~isempty(Result) && ~isempty(ct.cube)
       
        if isfield(Result,'RBE')
            Result.RBETruncated10Perc = Result.RBE;
            Result.RBETruncated10Perc(Result.physicalDose<0.1*...
                max(Result.physicalDose(:))) = 0;
        end

        DispInfo =fieldnames(Result);
        for i=1:size(DispInfo,1)
            
            if isstruct(Result.(DispInfo{i,1})) || isvector(Result.(DispInfo{i,1}))
                 Result = rmfield(Result,DispInfo{i,1});
                 DispInfo{i,2}=false;
            else
                %second dimension indicates if it should be plotted later on
                DispInfo{i,2}=true;
                DisablePlot = {'RBETruncated10Perc'};
                if strcmp(DispInfo{i,1},DisablePlot)
                    DispInfo{i,2}=false;
                end
                
                % determine units
                switch DispInfo{i,1}
                    case 'physicalDose'
                        DispInfo{i,3}  = '[Gy]';
                    case 'alpha'
                        DispInfo{i,3}  = '[Gy^{-1}]';
                    case 'beta'
                        DispInfo{i,3}  = '[Gy^{-2}]';
                    case 'RBExD'
                        DispInfo       = '[Gy(RBE)]';
                    otherwise
                        DispInfo{i,3}  = '[a.u.]';
                end
                
            end
        end

    set(handles.popupDisplayOption,'String',fieldnames(Result));
    set(handles.popupDisplayOption,'Value',find(strcmp(handles.SelectedDisplayOption,fieldnames(Result))));

    end
end

%% set and get required variables

plane = get(handles.popupPlane,'Value');
slice = round(get(handles.sliderSlice,'Value'));
CutOffLevel = 0.03;

%% plot ct
 if ~isempty(ct) && get(handles.popupTypeOfPlot,'Value')==1
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
if handles.State >2 &&  get(handles.popupTypeOfPlot,'Value')== 1

        mVolume = getfield(Result,handles.SelectedDisplayOption);
        % make sure to exploit full color range 
        mVolume(Result.physicalDose<...
            CutOffLevel*max(Result.physicalDose(:)))=0;

    %     %% dose colorwash
        if ~isempty(mVolume)&& ~isvector(mVolume)

            if handles.maxDoseVal == 0
                handles.maxDoseVal = max(mVolume(:));
                set(handles.txtMaxDoseVal,'String',num2str(handles.maxDoseVal))
            end
            dose_rgb = mVolume./handles.maxDoseVal;
            dose_rgb(dose_rgb>1) = 1;
            % Save RGB indices for dose in zsliceÂ´s voxels.
            if plane == 1  % Coronal plane
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(slice,:,:))),jet);
            elseif plane == 2 % sagittal plane
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,slice,:))),jet);
            elseif plane == 3 % Axial plane
                dose_rgb = ind2rgb(uint8(63*squeeze(dose_rgb(:,:,slice))),jet);
            end
            % plot dose distribution    
            doseImageHandle = image('CData',dose_rgb,'Parent',handles.axesFig);
 
            if ~isempty(ct.cube)
                if plane == 1  % Coronal plane
                    if get(handles.radiobtnDose,'Value') 
                        set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.(handles.SelectedDisplayOption)(slice,:,:))...
                            >CutOffLevel*handles.maxDoseVal)) ;
                    end   
                    
                elseif plane == 2 % sagittal plane
                    if get(handles.radiobtnDose,'Value')
                        set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.(handles.SelectedDisplayOption)(:,slice,:))>...
                        CutOffLevel*handles.maxDoseVal));
                    end
                elseif plane == 3 % Axial plane
                    if get(handles.radiobtnDose,'Value')
                        if strcmp(get(handles.popupDisplayOption,'String'),'RBETruncated10Perc')
                            set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.(handles.SelectedDisplayOption)(:,:,slice))...
                                >0.1*handles.maxDoseVal)) ;
                        else
                            set(doseImageHandle,'AlphaData',  .6*double(squeeze(Result.(handles.SelectedDisplayOption)(:,:,slice))...
                                >CutOffLevel*handles.maxDoseVal));
                        end
                    else
                        set(doseImageHandle,'AlphaData',  0*double(squeeze(Result.Dose(:,:,slice))>...
                            CutOffLevel*max(Result.(handles.SelectedDisplayOption)(:)))) ;
                    end
                end

            end

            % plot colorbar
            v=version;
            if str2double(v(1:3))>=8.5
                cBarHandel = colorbar(handles.axesFig,'colormap',jet,'FontSize',defaultFontSize,'yAxisLocation','right');
            else
                cBarHandel = colorbar('peer',handles.axesFig,'FontSize',defaultFontSize,'yAxisLocation','right');
            end
            Idx = find(strcmp(handles.SelectedDisplayOption,DispInfo(:,1)));
            set(get(cBarHandel,'ylabel'),'String', [DispInfo{Idx,1} ' ' DispInfo{Idx,3} ],'fontsize',defaultFontSize);
            
            if isempty(strfind(handles.SelectedDisplayOption,'RBE'))
                set(cBarHandel,'YLim',[0 handles.maxDoseVal]);
                caxis(handles.axesFig,[0,handles.maxDoseVal])
            else
                set(cBarHandel,'YLim',[0 handles.maxDoseVal]);
                caxis(handles.axesFig,[0,handles.maxDoseVal])
            end



        end
        
    axes(handles.axesFig),hold on 

        %% plot iso dose lines
        if get(handles.radiobtnIsoDoseLines,'Value')
                colormap(jet)
                SpacingLower = 0.1;
                SpacingUpper = 0.05;
                vLow  = 0.1:SpacingLower:0.9;
                vHigh = 0.95:SpacingUpper:1.2;
                vLevels = [vLow vHigh];
               
                if handles.IsoDose.RefVal == 0
                   MaxVal = max(mVolume(:)); 
                   vLevels = (round((vLevels.*MaxVal)*100))/100;
                else
                   MaxVal = handles.IsoDose.RefVal; 
                   vLevels = round(((handles.IsoDose.Levels)/100).*MaxVal*100)/100;
                end
                
                if plane == 1  % Coronal plane
                    Slice=squeeze(mVolume(slice,:,:));
                    if sum(Slice(:))>1
                        [C,myContour] = contour(Slice,vLevels);
                    end
                elseif plane == 2 % Sagittal plane
                    Slice=squeeze(mVolume(:,slice,:));
                    if sum(Slice(:))>1
                        [C,myContour] = contour(Slice,vLevels);
                    end
                elseif plane == 3 % Axial plane
                    Slice=squeeze(mVolume(:,:,slice));
                    if sum(Slice(:))>1
                        hold on
                     [C,myContour] = contour(Slice,vLevels,'LevelListMode','manual','LineWidth',1.5);
                    end
                end

                 if sum(Slice(:))>1
                    caxis(handles.axesFig,[0, handles.maxDoseVal]);
                    clabel(C,myContour,vLevels,'LabelSpacing',150)
                    % turn off legend for this data set
                    hAnnotation = get(myContour,'Annotation');
                    hLegendEntry = get(hAnnotation','LegendInformation');
                    set(hLegendEntry,'IconDisplayStyle','off');     
                    if get(handles.radiobtnIsoDoseLinesLabels,'Value') == 0
                        set(myContour,'ShowText','off')
                    end
                 end
        end

end


%% plot VOIs

if get(handles.radiobtnContour,'Value') && get(handles.popupTypeOfPlot,'Value')==1 && handles.State>0
    colors = colorcube;
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
    set(myLegend,'FontSize',defaultFontSize);
    set(myLegend,'color','none');
    set(myLegend,'TextColor', [1 1 1]);
    legend boxoff
    warning('on','MATLAB:legend:PlotEmpty')
end

%% Set axis labels
if  plane == 3% Axial plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution.x:1000);
        set(handles.axesFig,'YTick',0:50/ct.resolution.y:1000);
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution.x);
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution.y);   
        xlabel('x [mm]','FontSize',defaultFontSize)
        ylabel('y [mm]','FontSize',defaultFontSize)
        title(['axial plane z = ' num2str(ct.resolution.z*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel('x [voxels]','FontSize',defaultFontSize)
        ylabel('y [voxels]','FontSize',defaultFontSize)
        title('axial plane','FontSize',defaultFontSize)
    end
elseif plane == 2 % Sagittal plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution.z:1000)
        set(handles.axesFig,'YTick',0:50/ct.resolution.y:1000)
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution.z)
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution.y)
        xlabel('z [mm]','FontSize',defaultFontSize);
        ylabel('y [mm]','FontSize',defaultFontSize);
        title(['sagittal plane x = ' num2str(ct.resolution.y*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel('z [voxels]','FontSize',defaultFontSize)
        ylabel('y [voxels]','FontSize',defaultFontSize)
        title('sagittal plane','FontSize',defaultFontSize);
    end
elseif plane == 1 % Coronal plane
    if ~isempty(pln)
        set(handles.axesFig,'XTick',0:50/ct.resolution.z:1000)
        set(handles.axesFig,'YTick',0:50/ct.resolution.x:1000)
        set(handles.axesFig,'XTickLabel',0:50:1000*ct.resolution.z)
        set(handles.axesFig,'YTickLabel',0:50:1000*ct.resolution.x)
        xlabel('z [mm]','FontSize',defaultFontSize)
        ylabel('x [mm]','FontSize',defaultFontSize)
        title(['coronal plane y = ' num2str(ct.resolution.x*slice) ' [mm]'],'FontSize',defaultFontSize)
    else
        xlabel('z [voxels]','FontSize',defaultFontSize)
        ylabel('x [voxels]','FontSize',defaultFontSize)
        title('coronal plane','FontSize',defaultFontSize)
    end
end


%% profile plot
if get(handles.popupTypeOfPlot,'Value')==2 && exist('Result','var')
    % set SAD
    fileName = [pln.radiationMode '_' pln.machine];
    try
        load(fileName);
        SAD = machine.meta.SAD;
    catch
        error(['Could not find the following machine file: ' fileName ]); 
    end
     
    % clear view and initialize some values
    cla(handles.axesFig,'reset')
    set(gca,'YDir','normal');
    ylabel('{\color{black}dose [Gy]}')
    cColor={'black','green','magenta','cyan','yellow','red','blue'};
    
    % Rotation around Z axis (table movement)
    inv_rotMx_XY_T = [ cosd(pln.gantryAngles(handles.SelectedBeam)) sind(pln.gantryAngles(handles.SelectedBeam)) 0;
                      -sind(pln.gantryAngles(handles.SelectedBeam)) cosd(pln.gantryAngles(handles.SelectedBeam)) 0;
                                                                  0 0 1];
    % Rotation around Y axis (Couch movement)
    inv_rotMx_XZ_T = [cosd(pln.couchAngles(handles.SelectedBeam)) 0 -sind(pln.couchAngles(handles.SelectedBeam));
                                                                0 1 0;
                      sind(pln.couchAngles(handles.SelectedBeam)) 0 cosd(pln.couchAngles(handles.SelectedBeam))];
    
    if strcmp(handles.ProfileType,'longitudinal')
        sourcePointBEV = [handles.profileOffset -SAD   0];
        targetPointBEV = [handles.profileOffset  SAD   0];
    elseif strcmp(handles.ProfileType,'lateral')
        sourcePointBEV = [-SAD handles.profileOffset   0];
        targetPointBEV = [ SAD handles.profileOffset   0];
    end
    
    rotSourcePointBEV = sourcePointBEV * inv_rotMx_XZ_T * inv_rotMx_XY_T;
    rotTargetPointBEV = targetPointBEV * inv_rotMx_XZ_T * inv_rotMx_XY_T;
    
    % perform raytracing on the central axis of the selected beam
    [~,l,rho,~,ix] = matRad_siddonRayTracer(pln.isoCenter,ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{ct.cube});
    d = [0 l .* rho{1}];
    % Calculate accumulated d sum.
    vX = cumsum(d(1:end-1));
    
    % this step is necessary if visualization is set to profile plot
    % and another optimization is carried out - set focus on GUI
    figHandles = get(0,'Children');
    idxHandle = [];
    if ~isempty(figHandles)
        v=version;
        if str2double(v(1:3))>= 8.5
            idxHandle = strcmp({figHandles(:).Name},'matRadGUI');
        else
            idxHandle = strcmp(get(figHandles,'Name'),'matRadGUI');
        end
    end
    figure(figHandles(idxHandle));
    
    % plot physical dose
    mPhysDose = Result.('physicalDose'); 
    PlotHandles{1} = plot(handles.axesFig,vX,mPhysDose(ix),'color',cColor{1,1},'LineWidth',3); hold on; 
    PlotHandles{1,2} ='physicalDose';
    ylabel(handles.axesFig,'dose in [Gy]');
    set(handles.axesFig,'FontSize',defaultFontSize);
    
    % plot counter
    Cnt=2;
    
    if isfield(Result,'RBE')
        
        %disbale specific plots
        %DispInfo{6,2}=0;
        %DispInfo{5,2}=0;
        %DispInfo{2,2}=0;
        
        % generate two lines for ylabel
        StringYLabel1 = '\fontsize{8}{\color{red}RBE x dose [Gy(RBE)] \color{black}dose [Gy] ';
        StringYLabel2 = '';
        for i=1:1:size(DispInfo,1)
            if DispInfo{i,2} 
                %physicalDose is already plotted and RBExD vs RBE is plotted later with plotyy
                if ~strcmp(DispInfo{i,1},'RBExDose') &&...
                   ~strcmp(DispInfo{i,1},'RBE') && ...
                   ~strcmp(DispInfo{i,1},'physicalDose')
               
                        mCube = Result.(DispInfo{i,1});
                        PlotHandles{Cnt,1} = plot(handles.axesFig,vX,mCube(ix),'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                        PlotHandles{Cnt,2} = DispInfo{i,1};
                        StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' DispInfo{i,1} ' ['  DispInfo{i,3} ']'];
                        Cnt = Cnt+1;
                end    
            end
        end
        StringYLabel2 = [StringYLabel2 '}'];
        % always plot RBExD against RBE
        mRBExDose = Result.('RBExDose');
        vBED = mRBExDose(ix);
        mRBE = Result.('RBE');
        vRBE = mRBE(ix);
        
        % plot biological dose against RBE
        [ax, PlotHandles{Cnt,1}, PlotHandles{Cnt+1,1}]=plotyy(handles.axesFig,vX,vBED,vX,vRBE,'plot');hold on;
        PlotHandles{Cnt,2}='RBExDose';
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
        if strcmp(cst{i,3},'TARGET') && tmpPrior >= cst{i,5}.Priority && tmpSize<numel(cst{i,4})
           linIdxTarget = unique(cst{i,4});
           tmpPrior=cst{i,5}.Priority;
           tmpSize=numel(cst{i,4});
           VOI = cst{i,2};
        end
    end
    
    str = sprintf('profile plot - central axis of %d beam gantry angle %d° couch angle %d°',...
        handles.SelectedBeam ,pln.gantryAngles(handles.SelectedBeam),pln.couchAngles(handles.SelectedBeam));
    h_title = title(handles.axesFig,str,'FontSize',defaultFontSize);
    pos = get(h_title,'Position');
    set(h_title,'Position',[pos(1)-40 pos(2) pos(3)])
    
    % plot target boundaries
    mTargetCube = zeros(size(ct.cube));
    mTargetCube(linIdxTarget) = 1;
    vProfile = mTargetCube(ix);
    WEPL_Target_Entry = vX(find(vProfile,1,'first'));
    WEPL_Target_Exit  = vX(find(vProfile,1,'last'));
    PlotHandles{Cnt,2} =[VOI ' boundary'];
    
    

    if ~isempty(WEPL_Target_Entry) && ~isempty(WEPL_Target_Exit)
        hold on
        PlotHandles{Cnt,1} = ...
        plot([WEPL_Target_Entry WEPL_Target_Entry],handles.axesFig.YLim,'--','Linewidth',3,'color','k');hold on
        plot([WEPL_Target_Exit WEPL_Target_Exit], handles.axesFig.YLim,'--','Linewidth',3,'color','k');hold on
      
    else
        PlotHandles{Cnt,1} =[];
    end
    
   Lines  = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),1);
   Labels = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),2);
   h=legend(handles.axesFig,[Lines{:}],Labels{:});
   set(h,'FontSize',defaultFontSize);
   xlabel('radiological depth [mm]','FontSize',8);  
   grid on, grid minor
   
end




% --- Executes on selection change in popupPlane.
function popupPlane_Callback(hObject, ~, handles)
% hObject    handle to popupPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupPlane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupPlane

% set slice slider
handles.plane = get(handles.popupPlane,'value');
try
    if handles.State > 0
        ct = evalin('base', 'ct');
        set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
                'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);
        if handles.State < 3
            set(handles.sliderSlice,'Value',round(size(ct.cube,handles.plane)/2));
        else
            pln = evalin('base','pln');
            
            if handles.plane == 1
                set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.x));
            elseif handles.plane == 2
                set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.y));
            elseif handles.plane == 3
                set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.z));
            end
            
        end
    end
catch
end
        
UpdatePlot(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupPlane_CreateFcn(hObject, ~, ~)
% hObject    handle to popupPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderSlice_Callback(~, ~, handles)
% hObject    handle to sliderSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
UpdatePlot(handles)

% --- Executes during object creation, after setting all properties.
function sliderSlice_CreateFcn(hObject, ~, ~)
% hObject    handle to sliderSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in radiobtnContour.
function radiobtnContour_Callback(~, ~, handles)
% hObject    handle to radiobtnContour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnContour
UpdatePlot(handles)

% --- Executes on button press in radiobtnDose.
function radiobtnDose_Callback(~, ~, handles)
% hObject    handle to radiobtnDose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnDose
UpdatePlot(handles)

% --- Executes on button press in radiobtnIsoDoseLines.
function radiobtnIsoDoseLines_Callback(~, ~, handles)
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

try
    % indicate that matRad is busy
    % change mouse pointer to hour glass 
    Figures = findobj('type','figure');
    set(Figures, 'pointer', 'watch'); 
    drawnow;
    % disable all active objects
    InterfaceObj = findobj(Figures,'Enable','on');
    set(InterfaceObj,'Enable','off');
    % wait until the table is updated
    pause(0.1);
    uiTable_CellEditCallback(hObject,[],handles);
    pause(0.3);

    % if a critical change to the cst has been made which affects the dij matrix
    if handles.DijCalcWarning == true

        choice = questdlg('Overlap priorites of OAR constraints have been edited, a new OAR VOI was added or a critical row constraint was deleted. A new Dij calculation might be necessary.', ...
        'Title','Cancel','Calculate Dij then Optimize','Optimze directly','Optimze directly');

        switch choice
            case 'Cancel'
                return
            case 'Calculate dij again and optimize'
                handles.DijCalcWarning = false;
                btnCalcDose_Callback(hObject, eventdata, handles)      
            case 'Optimze directly'
                handles.DijCalcWarning = false;       
        end
    end

    % get optimization parameters from GUI
    Param.numOfIter = str2double(get(handles.editNumIter,'String'));
    Param.prec      = str2double(get(handles.txtPrecisionOutput,'String'));
    BioOptType      = get(handles.btnTypBioOpt,'String');
    
    pln = evalin('base','pln');
    ct  = evalin('base','ct');
    % optimize
    resultGUI = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),pln,1,Param,BioOptType);
    assignin('base','resultGUI',resultGUI);

    %set some values
    if handles.plane == 1
        set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.x));
    elseif handles.plane == 2
        set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.y));
    elseif handles.plane == 3
        set(handles.sliderSlice,'Value',ceil(pln.isoCenter(1,handles.plane)/ct.resolution.z));
    end

    handles.State = 3;
    handles.SelectedDisplayOptionIdx = 1;
    handles.SelectedDisplayOption='physicalDose';
    handles.SelectedBeam = 1;
    UpdatePlot(handles);
    UpdateState(handles);

catch 
    handles = showError(handles,'OptimizeCallback: Could not optimize'); 
    guidata(hObject,handles);
    return;
end


% perform sequencing and DAO
try
    
    %% sequencing
    if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
    %   resultGUI = matRad_xiaLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
    %       ,get(handles.editSequencingLevel,'Value'));
        resultGUI = matRad_engelLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
            ,str2double(get(handles.editSequencingLevel,'String')));
        assignin('base','resultGUI',resultGUI);
    end
catch
   handles = showError(handles,'OptimizeCallback: Could not perform sequencing'); 
   guidata(hObject,handles);
   return;
end

try
    %% DAO
    if strcmp(pln.radiationMode,'photons') && pln.runDAO
       resultGUI = matRad_directApertureOptimization(evalin('base','dij'),evalin('base','cst')...
           ,resultGUI.apertureInfo,resultGUI,pln,1);
       assignin('base','resultGUI',resultGUI);
    end
    
    if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
        matRad_visApertureInfo(resultGUI.apertureInfo);
    end

catch
   handles = showError(handles,'OptimizeCallback: Could not perform direct aperture optimization'); 
   guidata(hObject,handles);
   return;
end

% change state from busy to normal
set(Figures, 'pointer', 'arrow');
set(InterfaceObj,'Enable','on');
    
guidata(hObject,handles);



% the function CheckValidityPln checks if the provided plan is valid so
% that it can be used further on for optimization
function FlagValid = CheckValidityPln(cst)

FlagValid = true;
%check if mean constraint is always used in combination
for i = 1:size(cst,1)
   if ~isempty(cst{i,6})
        if ~isempty(strfind([cst{i,6}.type],'mean')) && isempty(strfind([cst{i,6}.type],'square'))
             FlagValid = false;
             warndlg('mean constraint needs to be defined in addition to a second constraint (e.g. squared deviation)');
             break      
        end
   end
end


% --- Executes on selection change in popupTypeOfPlot
function popupTypeOfPlot_Callback(hObject, ~, handles)
% hObject    handle to popupTypeOfPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 % intensity plot
if get(hObject,'Value') == 1  
    
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
    
% profile plot
elseif get(hObject,'Value') == 2
    
    if handles.State > 0
        if length(str2double(strsplit(get(handles.editGantryAngle,'String'),' '))) > 1
            
            set(handles.sliderBeamSelection,'Enable','on');
            handles.SelectedBeam = 1;
            pln = evalin('base','pln');
            set(handles.sliderBeamSelection,'Min',handles.SelectedBeam,'Max',pln.numOfBeams,...
                'Value',handles.SelectedBeam,...
                'SliderStep',[1/(pln.numOfBeams-1) 1/(pln.numOfBeams-1)],...
                'Enable','on');

        else
            handles.SelectedBeam = 1;
        end
    
        handles.profileOffset = get(handles.sliderOffset,'Value');

        vMinMax = [-100 100];
        vRange = sum(abs(vMinMax));

        ct = evalin('base','ct');
        if strcmp(get(handles.btnProfileType,'String'),'lateral')
            SliderStep = vRange/ct.resolution.x;       
        else
            SliderStep = vRange/ct.resolution.y;  
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
function popupTypeOfPlot_CreateFcn(hObject, ~, ~)
% hObject    handle to popupTypeOfPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupDisplayOption.
function popupDisplayOption_Callback(hObject, ~, handles)
% hObject    handle to popupDisplayOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupDisplayOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDisplayOption
content = get(hObject,'String');
handles.SelectedDisplayOption = content{get(hObject,'Value'),1};
handles.SelectedDisplayOptionIdx = get(hObject,'Value');
handles.maxDoseVal = 0;
guidata(hObject, handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function popupDisplayOption_CreateFcn(hObject, ~, ~)
% hObject    handle to popupDisplayOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderBeamSelection_Callback(hObject, ~, handles)
% hObject    handle to sliderBeamSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.SelectedBeam = round(get(hObject,'Value'));
set(hObject, 'Value', handles.SelectedBeam);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function sliderBeamSelection_CreateFcn(hObject, ~, ~)
% hObject    handle to sliderBeamSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnProfileType.
function btnProfileType_Callback(hObject, ~, handles)
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

% displays the cst in the GUI
function setCstTable(handles,cst)


columnname = {'VOI','VOI Type','Priority','Obj Func','Penalty','Parameter'};
columnformat = {cst(:,2)',{'OAR','TARGET'},'numeric',...
       {'square underdosing','square overdosing','square deviation', 'mean', 'EUD'},...
       'numeric','numeric'};
   
numOfObjectives = 0;
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        numOfObjectives = numOfObjectives + numel(cst{i,6});
    end
end

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
           case 'mean'
               data{Counter,6}='';
       
           case {'square underdosing','square overdosing','square deviation','EUD'}
               %Dose
               data{Counter,6}=cst{i,6}(j).parameter(1,2); 
       end
   
       Counter = Counter +1;
       end
   end
   
end

set(handles.uiTable,'ColumnName',columnname);
set(handles.uiTable,'ColumnFormat',columnformat);
set(handles.uiTable,'ColumnEditable',[true true true true true true true]);
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
                if isempty(data{j,1}) || ~isempty(strfind(data{j,1}, 'Select'))
                    FlagValidParameters=false;
                else
                    NewCst{Cnt,1}=data{j,1}; 
                end
                %VOI Type
                if isempty(data{j,2})|| ~isempty(strfind(data{j,2}, 'Select'))
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
            if isempty(data{j,4}) ||~isempty(strfind(data{j,4}, 'Select'))
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
                    if isempty(data{j,6})
                       FlagValidParameters=false;
                    else
                        NewCst{Cnt,4}(CntObjF,1).parameter(1,2) = (data{j,6});
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
                          if ischar(data{j,6})
                              NewCst{Cnt,4}(CntObjF,1).parameter(1,2) = str2double(data{j,6});
                          else
                              NewCst{Cnt,4}(CntObjF,1).parameter(1,2) = double(data{j,6});
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


% --- Executes on button press in btnuiTableAdd.
function btnuiTableAdd_Callback(hObject, ~, handles)
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

%handles.State=1;
guidata(hObject,handles);
UpdateState(handles);


% --- Executes on button press in btnuiTableDel.
function btnuiTableDel_Callback(hObject, eventdata, handles)
% hObject    handle to btnuiTableDel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.uiTable, 'data');
Index = get(handles.uiTable,'UserData');
mask = (1:size(data,1))';
mask(Index(:,1))=[];
cst = evalin('base','cst');
% if all rows have been deleted or a target voi was removed
if size(data,1)==1 || strcmp(data(Index(1),2),'TARGET')
    handles.State=1;
end

try
    Idx = find(strcmp(cst(:,2),data(Index(1),1)));
    % if OAR was removed then show a warning 
    if strcmp(data(Index(1),2),'OAR') && length(cst{Idx,6})<=1
      handles.DijCalcWarning =true;
    end
catch 
end
data=data(mask,:);
set(handles.uiTable,'data',data);
guidata(hObject,handles);
UpdateState(handles);
btnTableSave_Callback(hObject, eventdata, handles);

% --- Executes when selected cell(s) is changed in uiTable.
function uiTable_CellSelectionCallback(hObject, eventdata, ~)
% hObject    handle to uiTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
index = eventdata.Indices;
    if any(index)             %loop necessary to surpress unimportant errors.
        set(hObject,'UserData',index);      
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


% get table data and current index of cell
if isempty(eventdata)
    data =get(handles.uiTable,'Data');
    Index = get(handles.uiTable,'UserData');

    if ~isempty(Index) && size(Index,1)==1
        % if this callback was invoked by calculate dij button, eventdata is empty
        % and needs to be set manually
        try 
            % if row gots deleted then index is pointing to non existing
            % data
            if size(data,1)<Index(1,1)
                Index(1,1)=1;
                Index(1,2)=1;
            end
            eventdata.Indices(1) = Index(:,1);
            eventdata.Indices(2) = Index(:,2);
            eventdata.PreviousData = 1;
            eventdata.NewData = data{Index(1),Index(2)};
        catch
        end
    else
        return
    end
else
    data = get(hObject,'Data');
    
    
end

%% if VOI, VOI Type or Overlap was changed
if ~strcmp(eventdata.NewData,eventdata.PreviousData)
    if eventdata.Indices(2) == 1 || eventdata.Indices(2) == 2 ...
            || eventdata.Indices(2) == 3
        
        handles.TableChanged = true;
          %% if a new target is defined set state to one
         if eventdata.Indices(2) == 2 && strcmp('TARGET',eventdata.NewData)
              handles.State=1;
         end
        
        %% if overlap priority of target changed
         if eventdata.Indices(2) == 3 && strcmp('TARGET',data(eventdata.Indices(1),2))
              handles.State=1;
         end
         
        %% if overlap priority of OAR changed
         if eventdata.Indices(2) == 3 && strcmp('OAR',data(eventdata.Indices(1),2))
              handles.DijCalcWarning = true;
         end
         
         %% check if new OAR was added
         cst=evalin('base','cst');
         Idx = ~cellfun('isempty',cst(:,6));
         
         if sum(strcmp(cst(Idx,2),eventdata.NewData))==0 
             handles.DijCalcWarning = true;
         end
        
        
    else
        if handles.State ==3 && handles.TableChanged == false
            handles.State=2;
        end
    end
end
%% if VOI Type was changed -> check if objective function still makes sense
if eventdata.Indices(2) == 2
   
    if strcmp(eventdata.NewData,'OAR')
        
        if sum(strcmp({'square deviation','square underdosing'},data{eventdata.Indices(1),4}))>0
            data{eventdata.Indices(1),4} = 'square overdosing';
        end
        
    else
        
        if sum(strcmp({'EUD','mean'},data{eventdata.Indices(1),4}))>0
            data{eventdata.Indices(1),4} = 'square deviation';
        end
        
    end
end
%% if objective function was changed -> check if VOI Type still makes sense
if eventdata.Indices(2) == 4
   
    if strcmp(eventdata.NewData,'square deviation') || strcmp(eventdata.NewData,'square underdosing') 
        
        if strcmp('OAR',data{eventdata.Indices(1),2})
            data{eventdata.Indices(1),4} = 'square overdosing';
        end
        
    elseif strcmp(eventdata.NewData,'EUD') || strcmp(eventdata.NewData,'mean')
        
        if strcmp('TARGET',data{eventdata.Indices(1),2})
            data{eventdata.Indices(1),4} = 'square deviation';
        end
        
    end
end


%% check if input is a valid
%check if overlap,penalty and and parameters are numbers
if eventdata.Indices(2) == 3  || eventdata.Indices(2) == 5 || eventdata.Indices(2) == 6 ...
        && ~isempty(eventdata.NewData)
    if CheckValidity(eventdata.NewData) == false
            data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end
end

%% if objective function is set to mean --> set dose cell to empty
if eventdata.Indices(2)==4 && strcmp(eventdata.NewData,'mean')
     data{eventdata.Indices(1),6}='';
elseif strcmp(data{eventdata.Indices(1),4},'mean') && eventdata.Indices(2) == 6
     data{eventdata.Indices(1),6}='';
end


%% if VOI was changed then change VOI type and overlap according to new VOI
if eventdata.Indices(2) == 1 && eventdata.Indices(1) == size(data,1)
    for i = 1:size(data,1)
        if strcmp(eventdata.NewData,data{i,1})
           data{eventdata.Indices(1),2}=data{i,2};
           data{eventdata.Indices(1),3}=data{i,3};
        end
    end
    
end

%% set VOI type and priority according to existing definitions
for i=1:size(data,1)
    if i~=eventdata.Indices(1) && strcmp(data(i,1),data(eventdata.Indices(1)))
        data{i,2} = data{eventdata.Indices(1),2};
        data{i,3} = data{eventdata.Indices(1),3};
    end
end


if isnan(eventdata.NewData)
    data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
end


set(handles.txtInfo,'String','plan changed');
set(handles.uiTable,'data',data);
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

% enables/ disables buttons according to the current state      
function UpdateState(handles)

 switch handles.State
     
     case 0
      
      set(handles.txtInfo,'String','no data loaded');
      set(handles.btnCalcDose,'Enable','off');
      set(handles.btnOptimize ,'Enable','off');
      set(handles.btnDVH,'Enable','off');
      set(handles.btnSequencing,'Enable','off');
      
     case 1
      pln = evalin('base','pln');   
      set(handles.txtInfo,'String','ready for dose calculation');
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','off');  
      set(handles.btnDVH,'Enable','off');
      if strcmp(pln.radiationMode,'photons')
          set(handles.btnSequencing,'Enable','off');
      end
     case 2
      pln = evalin('base','pln');
      set(handles.txtInfo,'String','ready for optimization');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.btnDVH,'Enable','off');
      if strcmp(pln.radiationMode,'photons')
          set(handles.btnSequencing,'Enable','off');
      end
     
     case 3

      set(handles.txtInfo,'String','plan is optimized');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.btnDVH,'Enable','on');
      
      pln = evalin('base','pln');
      
      if strcmp(pln.radiationMode,'photons')
          set(handles.btnSequencing,'Enable','on');
      end
 end

 
 
% fill GUI elements with plan information
function setPln(handles)
pln=evalin('base','pln');
set(handles.editBixelWidth,'String',num2str(pln.bixelWidth));
set(handles.editFraction,'String',num2str(pln.numOfFractions));

if isfield(pln,'isoCenter')
    set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.isoCenter*10))./10), '\s+', ' '));
end
set(handles.editGantryAngle,'String',num2str((pln.gantryAngles)));
set(handles.editCouchAngle,'String',num2str((pln.couchAngles)));
set(handles.popupRadMode,'Value',find(strcmp(get(handles.popupRadMode,'String'),pln.radiationMode)));
set(handles.popUpMachine,'Value',find(strcmp(get(handles.popUpMachine,'String'),pln.machine)));

if strcmp(pln.bioOptimization,'effect') || strcmp(pln.bioOptimization,'RBExD') ... 
                            && strcmp(pln.radiationMode,'carbon')  
    set(handles.radbtnBioOpt,'Value',1);
    set(handles.radbtnBioOpt,'Enable','on');
    set(handles.btnTypBioOpt,'Enable','on');
    set(handles.btnTypBioOpt,'String',pln.bioOptimization);
else
    set(handles.radbtnBioOpt,'Value',0);
    set(handles.radbtnBioOpt,'Enable','off');
    set(handles.btnTypBioOpt,'Enable','off');
end
%% enable sequencing and DAO button if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons') && pln.runSequencing
    set(handles.btnRunSequencing,'Enable','on');
    set(handles.btnRunSequencing,'Value',1);
elseif strcmp(pln.radiationMode,'photons') && ~pln.runSequencing
    set(handles.btnRunSequencing,'Enable','on');
    set(handles.btnRunSequencing,'Value',0);
else
    set(handles.btnRunSequencing,'Enable','off');
end
%% enable DAO button if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons') && pln.runDAO
    set(handles.btnRunDAO,'Enable','on');
    set(handles.btnRunDAO,'Value',1);
elseif strcmp(pln.radiationMode,'photons') && ~pln.runDAO
    set(handles.btnRunDAO,'Enable','on');
    set(handles.btnRunDAO,'Value',0);
else
    set(handles.btnRunDAO,'Enable','off');
end
%% enable stratification level input if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons')
    set(handles.txtSequencing,'Enable','on');
    set(handles.editSequencingLevel,'Enable','on');
else
    set(handles.txtSequencing,'Enable','off');
    set(handles.editSequencingLevel,'Enable','off');
end

 
% get pln file form GUI     
function getPln(handles)

pln.bixelWidth      = parseStringAsNum(get(handles.editBixelWidth,'String'),false); % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = parseStringAsNum(get(handles.editGantryAngle,'String'),true); % [°]
pln.couchAngles     = parseStringAsNum(get(handles.editCouchAngle,'String'),true); % [°]
pln.numOfBeams      = numel(pln.gantryAngles);
try
    ct=evalin('base','ct');
    pln.numOfVoxels     = numel(ct.cube);
    pln.voxelDimensions = size(ct.cube);
catch
end
pln.numOfFractions  = parseStringAsNum(get(handles.editFraction,'String'),false);
contents            = get(handles.popupRadMode,'String'); 
pln.radiationMode   = contents{get(handles.popupRadMode,'Value')}; % either photons / protons / carbon
contents            = get(handles.popUpMachine,'String'); 
pln.machine         = contents{get(handles.popUpMachine,'Value')}; 

if (logical(get(handles.radbtnBioOpt,'Value')) && strcmp(pln.radiationMode,'carbon'))
    pln.bioOptimization =get(handles.btnTypBioOpt,'String');
else
     pln.bioOptimization = 'none';
end

pln.runSequencing = logical(get(handles.btnRunSequencing,'Value'));
pln.runDAO = logical(get(handles.btnRunDAO,'Value'));

try
    cst= evalin('base','cst');
    if sum(strcmp('TARGET',cst(:,3)))>0 && get(handles.checkIsoCenter,'Value')
       pln.isoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct')); 
    else
       pln.isoCenter = str2double(get(handles.editIsoCenter,'String'));
    end
catch
    warning('couldnt set isocenter in getPln function')
end

handles.pln = pln;
assignin('base','pln',pln);

function Number = parseStringAsNum(stringIn,isVector)
 Number = str2num(stringIn);
 if isempty(Number) || length(Number)>1 && ~isVector
     warndlg('could not parse all plan parameters'); 
end

% --- Executes on button press in btnTableSave.
function btnTableSave_Callback(~, ~, handles)
% hObject    handle to btnTableSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getCstTable(handles);
if get(handles.checkIsoCenter,'Value')
    pln = evalin('base','pln'); 
    pln.isoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct')); 
    set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.isoCenter*10))./10), '\s+', ' '));
    assignin('base','pln',pln);
end
getPln(handles);


% --- Executes on button press in btnDVH.
function btnDVH_Callback(~, ~, ~)
% hObject    handle to btnDVH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
matRad_calcDVH(evalin('base','resultGUI'),evalin('base','cst'))


% --- Executes on selection change in listBoxCmd.
function listBoxCmd_Callback(hObject, ~, ~)
% hObject    handle to listBoxCmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listBoxCmd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listBoxCmd
numLines = size(get(hObject,'String'),1);
set(hObject, 'ListboxTop', numLines);

% --- Executes during object creation, after setting all properties.
function listBoxCmd_CreateFcn(hObject, ~, ~)
% hObject    handle to listBoxCmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobtnIsoDoseLinesLabels.
function radiobtnIsoDoseLinesLabels_Callback(~, ~, handles)
% hObject    handle to radiobtnIsoDoseLinesLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnIsoDoseLinesLabels
UpdatePlot(handles);


% --- Executes on slider movement.
function sliderOffset_Callback(hObject, ~, handles)
% hObject    handle to sliderOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.profileOffset = get(hObject,'Value');
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function sliderOffset_CreateFcn(hObject, ~, ~)
% hObject    handle to sliderOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function About_Callback(~, ~, ~)
% hObject    handle to About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'https://github.com/e0404/matRad/' 'email: matrad@dkfz.de'},'About');




% --- Executes on button press in btnRefresh.
function btnRefresh_Callback(hObject, ~, handles)
% hObject    handle to btnRefresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%parse variables from base workspace
AllVarNames = evalin('base','who');
handles.State = 0;

if ~isempty(AllVarNames)
    try
         if  sum(ismember(AllVarNames,'ct')) > 0
            % do nothing
        else
             handles = showError(handles,'BtnRefreshCallback: cst struct is missing');
         end

        if  sum(ismember(AllVarNames,'cst')) > 0
            setCstTable(handles,evalin('base','cst'));
        else
             handles = showError(handles,'BtnRefreshCallback: cst struct is missing');
        end

        if sum(ismember(AllVarNames,'pln')) > 0
             setPln(handles);
        else
            getPln(handles);
        end
       handles.State = 1;

    catch
        handles = showError(handles,'BtnRefreshCallback: Could not load ct/cst/pln');
        guidata(hObject,handles);
        return;
    end

    % check if a optimized plan was loaded
    if sum(ismember(AllVarNames,'stf')) > 0  && sum(ismember(AllVarNames,'dij')) > 0
        handles.State = 2;
    end

    if sum(ismember(AllVarNames,'resultGUI')) > 0
        handles.State = 3;
        handles.State = 3;
        handles.SelectedDisplayOptionIdx = 1;
        handles.SelectedDisplayOption='physicalDose';
        handles.SelectedBeam = 1;
    end

    if handles.State > 0
        ct = evalin('base','ct');
        set(handles.sliderSlice,'Min',1,'Max',size(ct.cube,handles.plane),...
                'Value',ceil(size(ct.cube,handles.plane)/2),...
                'SliderStep',[1/(size(ct.cube,handles.plane)-1) 1/(size(ct.cube,handles.plane)-1)]);      
    end

end
guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);





% --------------------------------------------------------------------
function toolbarSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbarSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save Plan and Table to workspace
btnTableSave_Callback(hObject, eventdata, handles);

try

    if handles.State > 0
        ct  = evalin('base','ct');
        cst = evalin('base','cst');
        pln = evalin('base','pln');
    end

    if handles.State > 1
       stf = evalin('base','stf');
       dij = evalin('base','dij');
    end

    if handles.State > 2
       resultGUI = evalin('base','resultGUI');
    end

    switch handles.State
        case 1
            uisave({'cst','ct','pln'});
        case 2
            uisave({'cst','ct','pln','stf','dij'});
        case 3
            uisave({'cst','ct','pln','stf','dij','resultGUI'});
    end

catch
     handles = showWarning(handles,'Could not save files'); 
end
guidata(hObject,handles);



function editNumIter_Callback(hObject, ~, ~)
% hObject    handle to editNumIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumIter as text
%        str2double(get(hObject,'String')) returns contents of editNumIter as a double
val=round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val));

% --- Executes during object creation, after setting all properties.
function editNumIter_CreateFcn(hObject, ~, ~)
% hObject    handle to editNumIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slicerPrecision_Callback(hObject, ~, handles)
% hObject    handle to slicerPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
set(hObject,'Value',round(val))
set(handles.txtPrecisionOutput,'String',['1e-' num2str(round(val))]);

% --- Executes during object creation, after setting all properties.
function slicerPrecision_CreateFcn(hObject, ~, ~)
% hObject    handle to slicerPrecision (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnTypBioOpt.
function btnTypBioOpt_Callback(hObject, ~, handles)
% hObject    handle to btnTypBioOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getPln(handles);
if strcmp(get(hObject,'String'),'effect')
    set(hObject,'String','RBExD');
else
    set(hObject,'String','effect');
end


% --- Executes on button press in btnSequencing.
function btnSequencing_Callback(hObject, ~, handles)
% hObject    handle to btnSequencing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% indicate that matRad is busy
% change mouse pointer to hour glass 
Figures = findobj('type','figure');
set(Figures, 'pointer', 'watch'); 
drawnow;
% disable all active objects
InterfaceObj = findobj(Figures,'Enable','on');
set(InterfaceObj,'Enable','off');

StratificationLevel = str2double(get(handles.editSequencingLevel,'String'));
resultGUI   = evalin('base','resultGUI');
pln         = evalin('base','pln');
resultGUI.w = resultGUI.wUnsequenced;

% perform sequencing 
try
%% sequencing
if strcmp(pln.radiationMode,'photons') && (pln.runSequencing || pln.runDAO)
%   resultGUI = matRad_xiaLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
%       ,StratificationLevel,1);
    resultGUI = matRad_engelLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
        ,StratificationLevel,1);
    assignin('base','resultGUI',resultGUI);
end
catch
   handles = showError(handles,'BtnSequencingCallback: Could not perform sequencing');
   guidata(hObject,handles);
end

%% DAO
try
if strcmp(pln.radiationMode,'photons') && pln.runDAO
   resultGUI = matRad_directApertureOptimization(evalin('base','dij'),evalin('base','cst')...
       ,resultGUI.apertureInfo,resultGUI,1);
   matRad_visApertureInfo(resultGUI.apertureInfo);
   assignin('base','resultGUI',resultGUI);
end
catch 
   handles = showError(handles,'BtnSequencingCallback: Could not perform direct aperture optimization');
   guidata(hObject,handles);
end

% change state from busy to normal
set(Figures, 'pointer', 'arrow');
set(InterfaceObj,'Enable','on');

guidata(hObject,handles);
UpdatePlot(handles);
UpdateState(handles);

  



function editSequencingLevel_Callback(~, ~, ~)
% hObject    handle to editSequencingLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSequencingLevel as text
%        str2double(get(hObject,'String')) returns contents of editSequencingLevel as a double


% --- Executes during object creation, after setting all properties.
function editSequencingLevel_CreateFcn(hObject, ~, ~)
% hObject    handle to editSequencingLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIsoCenter_Callback(hObject, ~, handles)
% hObject    handle to editIsoCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if handles.State > 0
        pln=evalin('base','pln');
        if ~get(handles.checkIsoCenter,'Value')
           TmpIsoCenter = str2double(get(hObject,'String'));
           if length(TmpIsoCenter) == 3
               pln.isoCenter =  TmpIsoCenter;
           end
        end
    end
    assignin('base','pln',pln);
    handles.State = 1;
    UpdateState(handles);
catch 
    handles = showError(handles,'EditIsoCenterCallback: Could not set iso center');
    guidata(hObject,handles);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editIsoCenter_CreateFcn(hObject, ~, ~)
% hObject    handle to editIsoCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkIsoCenter.
function checkIsoCenter_Callback(hObject, ~, handles)
% hObject    handle to checkIsoCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkIsoCenter
try
    if get(hObject,'Value')
        pln = evalin('base','pln'); 
        pln.isoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct')); 
        set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.isoCenter*10))./10), '\s+', ' '));
        assignin('base','pln',pln);
    end
catch
end


% --- Executes on button press in btnRunSequencing.
function btnRunSequencing_Callback(~, ~, handles)
% hObject    handle to btnRunSequencing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getPln(handles);

% --- Executes on button press in btnRunDAO.
function btnRunDAO_Callback(~, ~, handles)
% hObject    handle to btnRunDAO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
getPln(handles);



function txtMaxDoseVal_Callback(hObject, ~, handles)
% hObject    handle to txtMaxDoseVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxDoseVal as text
%        str2double(get(hObject,'String')) returns contents of txtMaxDoseVal as a double
handles.maxDoseVal =  str2double(get(hObject,'String'));
guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function txtMaxDoseVal_CreateFcn(hObject, ~, ~)
% hObject    handle to txtMaxDoseVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in btnSetIsoDoseLevels.
function btnSetIsoDoseLevels_Callback(hObject, ~, handles)
% hObject    handle to btnSetIsoDoseLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Enter a reference dose. 100% correspond to a dose of [Gy]. Enter 0 to show original levels:','Provide percental iso dose levles (e.g. 95 105). Please enter space-separated numbers:'};
def = {'60','20 40 60 80 90 95 100 105 110'};
Input = inputdlg(prompt,'Set iso dose levels ', [1 50],def);
try
if ~isempty(Input)
     handles.IsoDose.RefVal = str2double(Input{1,:});
     handles.IsoDose.Levels = sort(str2double(Input{2,:})); 
     if length(handles.IsoDose.Levels) == 1
         handles.IsoDose.Levels = [handles.IsoDose.Levels handles.IsoDose.Levels];
     end
else
     handles.IsoDose.RefVal = 0;
end
catch
    handles.IsoDose.RefVal = 0;
end
UpdatePlot(handles);
guidata(hObject,handles);



% --- Executes on selection change in popUpMachine.
function popUpMachine_Callback(hObject, ~, handles)
% hObject    handle to popUpMachine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popUpMachine contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popUpMachine
contents = cellstr(get(hObject,'String'));
checkRadiationComposition(handles);
if handles.State>0
    pln = evalin('base','pln');
    if handles.State>0 && ~strcmp(contents(get(hObject,'Value')),pln.machine)
        handles.State=1;
        UpdateState(handles);
        guidata(hObject,handles);
    end
   getPln(handles);
end


function Valid = checkRadiationComposition(handles)
Valid = true;
contents = cellstr(get(handles.popUpMachine,'String'));
Machine = contents{get(handles.popUpMachine,'Value')};
contents = cellstr(get(handles.popupRadMode,'String'));
radMod = contents{get(handles.popupRadMode,'Value')};

FoundFile = dir([radMod '_' Machine '.mat']);
if isdeployed
   FoundFile = [FoundFile, dir([ctfroot filesep 'matRad' filesep radMod '_' Machine '.mat'])];
end
if isempty(FoundFile)
    warndlg(['No base data available for machine: ' Machine]);
    Valid = false;
end

% --- Executes during object creation, after setting all properties.
function popUpMachine_CreateFcn(hObject, ~, ~)
% hObject    handle to popUpMachine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = showError(handles,Message)

if isfield(handles,'ErrorDlg')
    close(handles.ErrorDlg);
end
 handles.ErrorDlg = errordlg(Message);

 function handles = showWarning(handles,Message)

if isfield(handles,'WarnDlg')
    close(handles.WarnDlg);
end
 handles.WarnDlg = warndlg(Message);


% --------------------------------------------------------------------
function toolbarLoad_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbarLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
btnLoadMat_Callback(hObject, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
selection = questdlg('Do you really want to close matRad?',...
                     'Close matRad',...
                     'Yes','No','Yes');
 switch selection,
   case 'Yes',
     delete(hObject);
   case 'No'
     return
end
