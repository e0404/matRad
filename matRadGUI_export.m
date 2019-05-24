function varargout = matRadGUI_export(varargin)
% matRad GUI
%
% call
%      MATRADGUI_EXPORT, by itself, creates a new MATRADGUI_EXPORT or raises the existing
%      singleton*.
%
%      H = MATRADGUI_EXPORT returns the handle to a new MATRADGUI_EXPORT or the handle to
%      the existing singleton*.
%
%      MATRADGUI_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRADGUI_EXPORT.M with the given input arguments.
%
%      MATRADGUI_EXPORT('Property','Value',...) creates a new MATRADGUI_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matRadGUI_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matRadGUI_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

if ~isdeployed
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(genpath(matRadRootDir));
end

[env, versionString] = matRad_getEnvironment();

% abort for octave
switch env
     case 'MATLAB'
         
     case 'OCTAVE'
         fprintf(['matRad GUI not available for ' env ' ' versionString ' \n']);
         return;
     otherwise
         fprintf(['not yet tested']);
 end
        
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
                   'gui_OpeningFcn', @matRadGUI_export_OpeningFcn, ...
                   'gui_OutputFcn',  @matRadGUI_export_OutputFcn, ...
                   'gui_LayoutFcn',  @matRadGUI_export_LayoutFcn, ...
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

function handles = resetGUI(hObject, handles, varargin)
% enable opengl software rendering to visualize linewidths properly
if ispc
  opengl software
elseif ismac
  % opengl is not supported
end

% Choose default command line output for matRadGUI_export
handles.output = hObject;
%show matrad logo
axes(handles.axesLogo)
[im, ~, alpha] = imread('matrad_logo.png');
f = image(im);
axis equal off
set(f, 'AlphaData', alpha);
% show dkfz logo
axes(handles.axesDKFZ)
[im, ~, alpha] = imread('DKFZ_Logo.png');
f = image(im);
axis equal off;
set(f, 'AlphaData', alpha);

% turn off the datacursormode (since it does not allow to set scrolling
% callbacks
handles.dcm_obj = datacursormode(handles.figure1);
set(handles.dcm_obj,'DisplayStyle','window');
if strcmpi(get(handles.dcm_obj,'Enable'),'on')
  set(handles.dcm_obj,'Enable','off');
end
%Add the callback for the datacursor display
set(handles.dcm_obj,'UpdateFcn',@dataCursorUpdateFunction);

% set callback for scroll wheel function
set(gcf,'WindowScrollWheelFcn',@matRadScrollWheelFcn);

% change color of toobar but only the first time GUI is started
if handles.initialGuiStart
  hToolbar = findall(hObject,'tag','uitoolbar1');
  jToolbar = get(get(hToolbar,'JavaContainer'),'ComponentPeer');
  jToolbar.setBorderPainted(false);
  color = java.awt.Color.gray;
  % Remove the toolbar border, to blend into figure contents
  jToolbar.setBackground(color);
  % Remove the separator line between toolbar and contents
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  jFrame = get(handle(hObject),'JavaFrame');
  jFrame.showTopSeparator(false);
  jtbc = jToolbar.getComponents;
  for idx=1:length(jtbc)
      jtbc(idx).setOpaque(false);
      jtbc(idx).setBackground(color);
      for childIdx = 1 : length(jtbc(idx).getComponents)
          jtbc(idx).getComponent(childIdx-1).setBackground(color);
      end
  end
end


set(handles.legendTable,'String',{'no data loaded'});
% clear  VOIPlotFlag
if isfield(handles,'VOIPlotFlag')
  handles = rmfield(handles,'VOIPlotFlag');
end

%seach for availabes machines
handles.Modalities = {'photons','protons','carbon'};
for i = 1:length(handles.Modalities)
  pattern = [handles.Modalities{1,i} '_*'];
  if isdeployed
      Files = dir([ctfroot filesep 'matRad' filesep pattern]);
  else
      Files = dir([fileparts(mfilename('fullpath')) filesep pattern]);        
  end
  for j = 1:length(Files)
      if ~isempty(Files)
          MachineName = Files(j).name(numel(handles.Modalities{1,i})+2:end-4);
          if isfield(handles,'Machines')
              if sum(strcmp(handles.Machines,MachineName)) == 0
                handles.Machines{size(handles.Machines,2)+1} = MachineName;
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
handles.CutOffLevel            = 0.01; % relative cut off level for dose vis
handles.IsoDose.NewIsoDoseFlag = false;
handles.TableChanged           = false;
handles.State                  = 0;
handles.doseOpacity            = 0.6;
handles.IsoDose.Levels         = 0;
handles.dispWindow             = cell(3,2); % first dimension refers to the selected

% do not calculate / suggest isoCenter new by default
set(handles.checkIsoCenter, 'Value', 0);
set(handles.editIsoCenter,'Enable','on')

% suppose no ct cube in HU is available (because no ct could be available)
handles.cubeHUavailable = false;
% initial startup finished
handles.initialGuiStart = false;
guidata(hObject, handles);  
% eof resetGUI

function handles = reloadGUI(hObject, handles, ct, cst)
AllVarNames = handles.AllVarNames;

if ismember('ct',AllVarNames)
    % compute HU values
    if ~isfield(ct, 'cubeHU')
        ct = matRad_electronDensitiesToHU(ct);
        assignin('base','ct',ct);
    end
    if ~isfield(ct, 'cubeHU')
        handles.cubeHUavailable = false;
    else
        handles.cubeHUavailable = true;
    end
end

%set plan if available - if not create one
try 
    
    if ismember('pln',AllVarNames) && handles.State > 0
        % check if you are working with a valid pln
        pln = evalin('base','pln');
        if ~isfield(pln,'propStf')
            handles = showWarning(handles,'GUI OpeningFunc: Overwriting outdated pln format with default GUI pln');
            evalin('base','clear pln');
            getPlnFromGUI(handles);
        end
        setPln(handles);
    elseif handles.State > 0 
         getPlnFromGUI(handles);
         setPln(handles);
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

% precompute iso dose lines
if handles.State == 3
    handles = updateIsoDoseLineCache(handles);
end

%per default the first beam is selected for the profile plot
handles.selectedBeam = 1;
handles.plane = get(handles.popupPlane,'Value');
handles.DijCalcWarning = false;

planePermute = [2 1 3];

% set slice slider
if handles.State > 0
    if evalin('base','exist(''pln'',''var'')')
        currPln = evalin('base','pln');
        if sum(currPln.propStf.isoCenter(:)) ~= 0
            currSlice = ceil(currPln.propStf.isoCenter(1,planePermute(handles.plane))/ct.resolution.x);
        else 
            currSlice = ceil(ct.cubeDim(planePermute(handles.plane))/2);
        end
    else % no pln -> no isocenter -> use middle
        currSlice = ceil(ct.cubeDim(planePermute(handles.plane))/2);
    end
    set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(handles.plane),...
            'Value',currSlice,...
            'SliderStep',[1/(ct.cubeDim(handles.plane)-1) 1/(ct.cubeDim(handles.plane)-1)]);      
    
    % define context menu for structures
    for i = 1:size(cst,1)
        if cst{i,5}.Visible
            handles.VOIPlotFlag(i) = true;
        else
            handles.VOIPlotFlag(i) = false;
        end
    end
  else
    % reset slider when nothing is loaded
    set(handles.sliderSlice,'Min',0,'Max',1,'Value',0,'SliderStep',[1 1]);
end

%Initialize colormaps and windows
handles.doseColorMap = 'jet';
handles.ctColorMap   = 'bone';
handles.cMapSize     = 64;
handles.cBarChanged  = true;

%Set up the colormap selection box
availableColormaps = matRad_getColormap();
set(handles.popupmenu_chooseColormap,'String',availableColormaps);

currentCtMapIndex   = find(strcmp(availableColormaps,handles.ctColorMap));
currentDoseMapIndex = find(strcmp(availableColormaps,handles.doseColorMap));

if handles.State >= 1    
   set(handles.popupmenu_chooseColormap,'Value',currentCtMapIndex);
end

if handles.State >= 3
    set(handles.popupmenu_chooseColormap,'Value',currentDoseMapIndex);
end

% Update handles structure
handles.profileOffset = 0;
UpdateState(handles)

axes(handles.axesFig)

handles.rememberCurrAxes = false;
UpdatePlot(handles)
handles.rememberCurrAxes = true;
guidata(hObject, handles);  
% eof reloadGUI

% --- Executes just before matRadGUI_export is made visible.
function matRadGUI_export_OpeningFcn(hObject, ~, handles, varargin) 
%#ok<*DEFNU> 
%#ok<*AGROW>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRadGUI_export (see VARARGIN)

% variable to check whether GUI is opened or just refreshed / new data
% loaded, since resetGUI needs to distinguish at one point

handles.initialGuiStart = true;

%If devMode is true, error dialogs will include the full stack trace of the error
%If false, only the basic error message is shown (works for errors that
%handle the MException object)
handles.devMode = true;

set(handles.radiobtnPlan,'value',0);

handles = resetGUI(hObject, handles);

%% parse variables from base workspace
AllVarNames = evalin('base','who');
handles.AllVarNames = AllVarNames;
try
    if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
        ct  = evalin('base','ct');
        cst = evalin('base','cst');
        %cst = setCstTable(handles,cst);
        cst = generateCstTable(handles,cst);
        handles.State = 1;
        cst = matRad_computeVoiContoursWrapper(cst,ct);
        assignin('base','cst',cst);

    elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find cst file');
    elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find ct file');
    end
catch
   handles = showError(handles,'GUI OpeningFunc: Could not load ct and cst file');
end

if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
    handles = reloadGUI(hObject, handles, ct, cst);
else
    handles = reloadGUI(hObject, handles);
end

guidata(hObject, handles);


function Callback_StructVisibilty(source,~)

handles = guidata(findobj('Name','matRadGUI'));

contextUiChildren = get(get(handles.figure1,'UIContextMenu'),'Children');

Idx = find(strcmp(get(contextUiChildren,'Label'),get(source,'Label')));
if strcmp(get(source,'Checked'),'on')
    set(contextUiChildren(Idx),'Checked','off');
else
    set(contextUiChildren(Idx),'Checked','on');
end
%guidata(findobj('Name','matRadGUI_export'), handles);
UpdatePlot(handles);

Update
% --- Outputs from this function are returned to the command line.
function varargout = matRadGUI_export_OutputFcn(~, ~, handles) 
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

[FileName, FilePath] = uigetfile('*.mat');
if FileName == 0 % user pressed cancel --> do nothing.
    return;
end

handles = resetGUI(hObject, handles);

try 
    
    % delete existing workspace - parse variables from base workspace
    AllVarNames = evalin('base','who');
    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
    
    for i = 1:length(RefVarNames)  
        if sum(ismember(AllVarNames,RefVarNames{i}))>0
            evalin('base',['clear ', RefVarNames{i}]);
        end
    end

    % read new data
    load([FilePath FileName]);
    set(handles.legendTable,'String',{'no data loaded'});
    set(handles.popupDisplayOption,'String','no option available');
    
catch ME
    handles = showError(handles,'LoadMatFileFnc: Could not load *.mat file',ME); 

    guidata(hObject,handles);
    UpdateState(handles);
    UpdatePlot(handles);
    return
end

try
    generateCstTable(handles,cst);
    handles.TableChanged = false;
    set(handles.popupTypeOfPlot,'Value',1);
    cst = matRad_computeVoiContoursWrapper(cst,ct);

    assignin('base','ct',ct);
    assignin('base','cst',cst);
    handles.State = 1;
catch ME
    handles = showError(handles,'LoadMatFileFnc: Could not load *.mat file',ME); 
end

% check if a optimized plan was loaded
if exist('stf','var')
    assignin('base','stf',stf);
end
if exist('pln','var')
    assignin('base','pln',pln);
end
if exist('dij','var')
    assignin('base','dij',dij);
end
% if exist('stf','var') && exist('dij','var')
%     handles.State = 2;
% end

if exist('resultGUI','var')
    assignin('base','resultGUI',resultGUI);
    % handles.State = 3;
    % handles.SelectedDisplayOption ='physicalDose';
end

% recheck current workspace variables
AllVarNames = evalin('base','who');
handles.AllVarNames = AllVarNames;

if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
    handles = reloadGUI(hObject, handles, ct, cst);
else
    handles = reloadGUI(hObject, handles);
end

guidata(hObject,handles);

% --- Executes on button press in btnLoadDicom.
function btnLoadDicom_Callback(hObject, ~, handles)
% hObject    handle to btnLoadDicom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    % delete existing workspace - parse variables from base workspace
    set(handles.popupDisplayOption,'String','no option available');
    AllVarNames = evalin('base','who');
    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
    for i = 1:length(RefVarNames)  
        if sum(ismember(AllVarNames,RefVarNames{i}))>0
            evalin('base',['clear ', RefVarNames{i}]);
        end
    end
    handles.State = 0;
    matRad_importDicomGUI;
 
catch
   handles = showError(handles,'DicomImport: Could not import data'); 
end
UpdateState(handles);
guidata(hObject,handles);


% --- Executes on button press in btn_export.
function btn_export_Callback(hObject, eventdata, handles)
% hObject    handle to btn_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    matRad_exportGUI;
catch
    handles = showError(handles,'Could not export data'); 
end
UpdateState(handles);
guidata(hObject,handles);

function editBixelWidth_Callback(hObject, ~, handles)
% hObject    handle to editBixelWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBixelWidth as text
%        str2double(get(hObject,'String')) returns contents of editBixelWidth as a double
getPlnFromGUI(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

function editGantryAngle_Callback(hObject, ~, handles)
% hObject    handle to editGantryAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGantryAngle as text
%        str2double(get(hObject,'String')) returns contents of editGantryAngle as a double
getPlnFromGUI(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    UpdatePlot(handles);
    guidata(hObject,handles);
end

function editCouchAngle_Callback(hObject, ~, handles)
% hObject    handle to editCouchAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCouchAngle as text
%        str2double(get(hObject,'String')) returns contents of editCouchAngle as a double
getPlnFromGUI(handles);
if handles.State > 0
    handles.State = 1;
    UpdateState(handles);
    guidata(hObject,handles);
end

% --- Executes on selection change in popupRadMode.
function popupRadMode_Callback(hObject, eventdata, handles)
% hObject    handle to popupRadMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkRadiationComposition(handles);
contents      = cellstr(get(hObject,'String')); 
RadIdentifier = contents{get(hObject,'Value')};
contentPopUp  = get(handles.popMenuBioOpt,'String');
switch RadIdentifier
    case 'photons'

        set(handles.popMenuBioOpt,'Enable','off');
        ix = find(strcmp(contentPopUp,'none'));
        set(handles.popMenuBioOpt,'Value',ix);
        set(handles.btnSetTissue,'Enable','off');
        
        set(handles.btnRunSequencing,'Enable','on');
        set(handles.btnRunDAO,'Enable','on');
        set(handles.radiobutton3Dconf,'Enable','on');
        set(handles.txtSequencing,'Enable','on');
        set(handles.editSequencingLevel,'Enable','on');
        
    case 'protons'
        
        set(handles.popMenuBioOpt,'Enable','on');
        ix = find(strcmp(contentPopUp,'const_RBExD'));
        set(handles.popMenuBioOpt,'Value',ix);
        set(handles.btnSetTissue,'Enable','on');
        
        set(handles.btnSetTissue,'Enable','off');
        set(handles.btnRunSequencing,'Enable','off');
        set(handles.btnRunDAO,'Enable','off');
        set(handles.radiobutton3Dconf,'Enable','off');
        set(handles.txtSequencing,'Enable','off');
        set(handles.editSequencingLevel,'Enable','off');
        
    case 'carbon'

        set(handles.popMenuBioOpt,'Enable','on');
        ix = find(strcmp(contentPopUp,'LEMIV_RBExD'));
        set(handles.popMenuBioOpt,'Value',ix);
        set(handles.btnSetTissue,'Enable','on');
        
        set(handles.btnRunSequencing,'Enable','off');
        set(handles.btnRunDAO,'Enable','off');
        set(handles.radiobutton3Dconf,'Enable','off');
        set(handles.txtSequencing,'Enable','off');
        set(handles.editSequencingLevel,'Enable','off');
end

if handles.State > 0
    pln = evalin('base','pln');
    if handles.State > 0 && ~strcmp(contents(get(hObject,'Value')),pln.radiationMode)
        
        % new radiation modality is photons -> just keep physicalDose
        if strcmp(contents(get(hObject,'Value')),'photons')
            try  
            AllVarNames = evalin('base','who');
               if  ismember('resultGUI',AllVarNames)
                   resultGUI = evalin('base','resultGUI');
                   if isfield(resultGUI,'alpha');    resultGUI = rmfield(resultGUI,'alpha');   end
                   if isfield(resultGUI,'beta');     resultGUI = rmfield(resultGUI,'beta');    end
                   if isfield(resultGUI,'RBExDose'); resultGUI = rmfield(resultGUI,'RBExDose');end
                   if isfield(resultGUI,'RBE');      resultGUI = rmfield(resultGUI,'RBE');     end
                   assignin('base','resultGUI',resultGUI);
                   handles = updateIsoDoseLineCache(handles);
               end   
            catch
            end
        elseif strcmp(contents(get(hObject,'Value')),'protons')
            try  
            AllVarNames = evalin('base','who');
               if  ismember('resultGUI',AllVarNames)
                   resultGUI = evalin('base','resultGUI');
                   if isfield(resultGUI,'alpha'); resultGUI = rmfield(resultGUI,'alpha');end
                   if isfield(resultGUI,'beta');  resultGUI = rmfield(resultGUI,'beta'); end
                   if isfield(resultGUI,'RBE');   resultGUI = rmfield(resultGUI,'RBE');  end
                   assignin('base','resultGUI',resultGUI);
                   handles = updateIsoDoseLineCache(handles);
               end
            catch
            end
        end
      
        guidata(hObject,handles);
        UpdatePlot(handles);
       
        getPlnFromGUI(handles);
        handles.State = 1;
        UpdateState(handles);

    end
         
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
    Figures = gcf;%findobj('type','figure');
    set(Figures, 'pointer', 'watch'); 
    drawnow;
    % disable all active objects
    InterfaceObj = findobj(Figures,'Enable','on');
    set(InterfaceObj,'Enable','off');
    
    %pause(0.1);
    %uiTable_CellEditCallback(hObject,[],handles);
    %pause(0.3);

    %% get cst from table
    %if ~getCstTable(handles)
    %    return
    %end
    % read plan from gui and save it to workspace
    % gets also IsoCenter from GUI if checkbox is not checked
    getPlnFromGUI(handles);

    % get default iso center as center of gravity of all targets if not
    % already defined
    pln = evalin('base','pln');

    if length(pln.propStf.gantryAngles) ~= length(pln.propStf.couchAngles) 
        handles = showWarning(handles,'number of gantryAngles != number of couchAngles'); 
    end
    %%
    if ~checkRadiationComposition(handles);
        fileName = [pln.radiationMode '_' pln.machine];
        handles = showError(handles,errordlg(['Could not find the following machine file: ' fileName ]));
        guidata(hObject,handles);
        return;
    end

    %% check if isocenter is already set
    if ~isfield(pln.propStf,'isoCenter')
        handles = showWarning(handles,'no iso center set - using center of gravity based on structures defined as TARGET');
        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
        assignin('base','pln',pln);
    elseif ~get(handles.checkIsoCenter,'Value')
        if ~strcmp(get(handles.editIsoCenter,'String'),'multiple isoCenter')
            pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*str2num(get(handles.editIsoCenter,'String'));
        end
    end

catch ME
    handles = showError(handles,'CalcDoseCallback: Error in preprocessing!',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

% generate steering file
try 
    currPln = evalin('base','pln');
    % if we run 3d conf opt -> hijack runDao to trigger computation of
    % connected bixels
    if strcmp(pln.radiationMode,'photons') && get(handles.radiobutton3Dconf,'Value')
       currpln.propOpt.runDAO = true; 
    end
    stf = matRad_generateStf(evalin('base','ct'),...
                                     evalin('base','cst'),...
                                     currPln);
    assignin('base','stf',stf);
catch ME
    handles = showError(handles,'CalcDoseCallback: Error in steering file generation!',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

% carry out dose calculation
try
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(evalin('base','ct'),stf,pln,evalin('base','cst'));
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(evalin('base','ct'),stf,pln,evalin('base','cst'));
    end

    % assign results to base worksapce
    assignin('base','dij',dij);
    handles.State = 2;
    handles.TableChanged = false;
    UpdateState(handles);
    UpdatePlot(handles);
    guidata(hObject,handles);
catch ME
    handles = showError(handles,'CalcDoseCallback: Error in dose calculation!',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

% change state from busy to normal
set(Figures, 'pointer', 'arrow');
set(InterfaceObj,'Enable','on');

guidata(hObject,handles);  



%% plots ct and distributions
function UpdatePlot(handles)

%profile on;

axes(handles.axesFig);

% this is necessary to prevent multiple callbacks of update plot drawing on
% top of each other in matlab <2014
drawnow;

defaultFontSize = 8;
currAxes            = axis(handles.axesFig);
AxesHandlesCT_Dose  = gobjects(0);
AxesHandlesVOI      = cell(0);
AxesHandlesIsoDose  = gobjects(0);

if handles.State == 0
    cla reset
    return
elseif handles.State > 0
     ct  = evalin('base','ct');
     cst = evalin('base','cst');
     pln = evalin('base','pln');
end

%% state 3 indicates that a optimization has been performed
 AllVarNames = evalin('base','who');
if  ismember('resultGUI',AllVarNames)
    Result = evalin('base','resultGUI');
end

if exist('Result','var')
    if ~isempty(Result) && ~isempty(ct.cubeHU) && ~isfield(handles,'DispInfo')
        
        DispInfo = fieldnames(Result);
        
        for i = 1:size(DispInfo,1)
            
            % delete weight vectors in Result struct for plotting
            if isstruct(Result.(DispInfo{i,1})) || isvector(Result.(DispInfo{i,1}))
                 Result = rmfield(Result,DispInfo{i,1});
                 DispInfo{i,2}=false;
            else
                %second dimension indicates if it should be plotted 
                DispInfo{i,2} = true;
                % determine units
                if strfind(DispInfo{i,1},'physicalDose')
                    DispInfo{i,3} = '[Gy]';
                elseif strfind(DispInfo{i,1},'alpha')
                    DispInfo{i,3} = '[Gy^{-1}]';
                elseif strfind(DispInfo{i,1},'beta')
                    DispInfo{i,3} = '[Gy^{-2}]';
                elseif strfind(DispInfo{i,1},'RBExD')
                    DispInfo{i,3} = '[Gy(RBE)]';
                elseif strfind(DispInfo{i,1},'LET')
                    DispInfo{i,3} = '[keV/um]';
                else
                    DispInfo{i,3} = '[a.u.]';
                end
                DispInfo{i,4} = [];    % optional for the future: color range for plotting
                DispInfo{i,5} = [];    % optional for the future: min max values
            end
        end

        set(handles.popupDisplayOption,'String',fieldnames(Result));
        if sum(strcmp(handles.SelectedDisplayOption,fieldnames(Result))) == 0
            handles.SelectedDisplayOption = DispInfo{find([DispInfo{:,2}],1,'first'),1};
        end
        set(handles.popupDisplayOption,'Value',find(strcmp(handles.SelectedDisplayOption,fieldnames(Result))));

    end
end

%% set and get required variables
plane = get(handles.popupPlane,'Value');
slice = round(get(handles.sliderSlice,'Value'));
hold(handles.axesFig,'on');
if get(handles.popupTypeOfPlot,'Value')==1
    set(handles.axesFig,'YDir','Reverse');
end

%% Remove colorbar?
plotColorbarSelection = get(handles.popupmenu_chooseColorData,'Value');

if get(handles.popupTypeOfPlot,'Value')==2 || plotColorbarSelection == 1
    if isfield(handles,'cBarHandel')
        delete(handles.cBarHandel);
    end
    %The following seems to be necessary as MATLAB messes up some stuff 
    %with the handle storage
    ch = findall(gcf,'tag','Colorbar');
    if ~isempty(ch)
        delete(ch);
    end
end

selectIx = get(handles.popupmenu_chooseColorData,'Value');  

cla(handles.axesFig); 
%% plot ct - if a ct cube is available and type of plot is set to 1 and not 2; 1 indicate cube plotting and 2 profile plotting
if ~isempty(ct) && get(handles.popupTypeOfPlot,'Value')==1
    
    if selectIx == 3
        ctIx = 2;
    else
        ctIx = selectIx;
    end
    
    if isfield(ct, 'cube')
        plotCtCube = ct.cube;
    else
        plotCtCube = ct.cubeHU;
    end
    
    ctMap = matRad_getColormap(handles.ctColorMap,handles.cMapSize);
    
    if isempty(handles.dispWindow{ctIx,2})
        handles.dispWindow{ctIx,2} = [min(reshape([ct.cubeHU{:}],[],1)) max(reshape([ct.cubeHU{:}],[],1))];
    end

    if get(handles.radiobtnCT,'Value')
        [AxesHandlesCT_Dose(end+1),~,handles.dispWindow{ctIx,1}] = matRad_plotCtSlice(handles.axesFig,plotCtCube,1,plane,slice,ctMap,handles.dispWindow{ctIx,1});
        
        % plot colorbar? If 1 the user asked for the CT
        if plotColorbarSelection == 2 && handles.cBarChanged
            %Plot the colorbar
            handles.cBarHandel = matRad_plotColorbar(handles.axesFig,ctMap,handles.dispWindow{ctIx,1},'fontsize',defaultFontSize);
            %adjust lables
            if isfield(ct,'cubeHU')
                set(get(handles.cBarHandel,'ylabel'),'String', 'Hounsfield Units','fontsize',defaultFontSize);
            else
                set(get(handles.cBarHandel,'ylabel'),'String', 'Electron Density','fontsize',defaultFontSize);
            end
            % do not interprete as tex syntax
            set(get(handles.cBarHandel,'ylabel'),'interpreter','none');
        end
    end
end

%% plot dose cube
if handles.State >= 1 &&  get(handles.popupTypeOfPlot,'Value')== 1  && exist('Result','var')
        doseMap = matRad_getColormap(handles.doseColorMap,handles.cMapSize);
        doseIx  = 3;
        % if the selected display option doesn't exist then simply display
        % the first cube of the Result struct
        if ~isfield(Result,handles.SelectedDisplayOption)
            CubeNames = fieldnames(Result);
            handles.SelectedDisplayOption = CubeNames{1,1};
        end
        
        dose = Result.(handles.SelectedDisplayOption);
        
        % dose colorwash
        if ~isempty(dose) && ~isvector(dose)
           
            if isempty(handles.dispWindow{doseIx,2})
                handles.dispWindow{doseIx,2} = [min(dose(:)) max(dose(:))];   % set min and max dose values
            end
            
            if get(handles.radiobtnDose,'Value')
                [doseHandle,~,handles.dispWindow{doseIx,1}] = matRad_plotDoseSlice(handles.axesFig,dose,plane,slice,handles.CutOffLevel,handles.doseOpacity,doseMap,handles.dispWindow{doseIx,1});
                AxesHandlesCT_Dose(end+1)         = doseHandle;
            end            
                    
            % plot colorbar?
            if plotColorbarSelection > 2 && handles.cBarChanged
                %Plot the colorbar
                handles.cBarHandel = matRad_plotColorbar(handles.axesFig,doseMap,handles.dispWindow{selectIx,1},'fontsize',defaultFontSize);
                %adjust lables
                Idx = find(strcmp(handles.SelectedDisplayOption,DispInfo(:,1)));
                set(get(handles.cBarHandel,'ylabel'),'String', [DispInfo{Idx,1} ' ' DispInfo{Idx,3} ],'fontsize',defaultFontSize);
                % do not interprete as tex syntax
                set(get(handles.cBarHandel,'ylabel'),'interpreter','none');
            end
        end
        

        %% plot iso dose lines
        if get(handles.radiobtnIsoDoseLines,'Value')           
            plotLabels = get(handles.radiobtnIsoDoseLinesLabels,'Value') == 1;
            
            %Sanity Check for Contours, which actually should have been 
            %computed before calling UpdatePlot
            if ~isfield(handles.IsoDose,'Contours')
                try
                    handles.IsoDose.Contours = matRad_computeIsoDoseContours(dose,handles.IsoDose.Levels); 
                catch
                    %If the computation didn't work, we set the field to
                    %empty, which will force matRad_plotIsoDoseLines to use
                    %matlabs contour function instead of repeating the
                    %failing computation every time
                    handles.IsoDose.Contours = [];
                    warning('Could not compute isodose lines! Will try slower contour function!');
                end
            end
            AxesHandlesIsoDose = matRad_plotIsoDoseLines(handles.axesFig,dose,handles.IsoDose.Contours,handles.IsoDose.Levels,plotLabels,plane,slice,doseMap,handles.dispWindow{doseIx,1},'LineWidth',1.5);
        end
end

selectIx = get(handles.popupmenu_chooseColorData,'Value');
set(handles.txtMinVal,'String',num2str(handles.dispWindow{selectIx,2}(1,1)));
set(handles.txtMaxVal,'String',num2str(handles.dispWindow{selectIx,2}(1,2)));

%% plot VOIs
if get(handles.radiobtnContour,'Value') && get(handles.popupTypeOfPlot,'Value')==1 && handles.State>0
    AxesHandlesVOI = [AxesHandlesVOI matRad_plotVoiContourSlice(handles.axesFig,cst,ct,1,handles.VOIPlotFlag,plane,slice,[],'LineWidth',2)];
end

%% Set axis labels and plot iso center
matRad_plotAxisLabels(handles.axesFig,ct,plane,slice,defaultFontSize);

if get(handles.radioBtnIsoCenter,'Value') == 1 && get(handles.popupTypeOfPlot,'Value') == 1 && ~isempty(pln)
    hIsoCenterCross = matRad_plotIsoCenterMarker(handles.axesFig,pln,ct,plane,slice);
end

if get(handles.radiobtnPlan,'value') == 1 && ~isempty(pln)
    matRad_plotProjectedGantryAngles(handles.axesFig,pln,ct,plane);
end

% the following line ensures the plotting order (optional)
% set(gca,'Children',[AxesHandlesCT_Dose hIsoCenterCross AxesHandlesIsoDose  AxesHandlesVOI ]);    
  
%set axis ratio
ratios = [1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
set(handles.axesFig,'DataAspectRatioMode','manual');
if plane == 1 
      res = [ratios(3) ratios(2)]./max([ratios(3) ratios(2)]);  
      set(handles.axesFig,'DataAspectRatio',[res 1])
elseif plane == 2 % sagittal plane
      res = [ratios(3) ratios(1)]./max([ratios(3) ratios(1)]);  
      set(handles.axesFig,'DataAspectRatio',[res 1]) 
elseif  plane == 3 % Axial plane
      res = [ratios(2) ratios(1)]./max([ratios(2) ratios(1)]);  
      set(handles.axesFig,'DataAspectRatio',[res 1])
end


%% profile plot
if get(handles.popupTypeOfPlot,'Value') == 2 && exist('Result','var')
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
    
    % Rotate the system into the beam. 
    % passive rotation & row vector multiplication & inverted rotation requires triple matrix transpose                  
    rotMat_system_T = transpose(matRad_getRotationMatrix(pln.propStf.gantryAngles(handles.selectedBeam),pln.propStf.couchAngles(handles.selectedBeam)));
    
    if strcmp(handles.ProfileType,'longitudinal')
        sourcePointBEV = [handles.profileOffset -SAD   0];
        targetPointBEV = [handles.profileOffset  SAD   0];
    elseif strcmp(handles.ProfileType,'lateral')
        sourcePointBEV = [-SAD handles.profileOffset   0];
        targetPointBEV = [ SAD handles.profileOffset   0];
    end
    
    rotSourcePointBEV = sourcePointBEV * rotMat_system_T;
    rotTargetPointBEV = targetPointBEV * rotMat_system_T;
    
    % perform raytracing on the central axis of the selected beam, use unit
    % electron density for plotting against the geometrical depth
    [~,l,rho,~,ix] = matRad_siddonRayTracer(pln.propStf.isoCenter(handles.selectedBeam,:),ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{0*ct.cubeHU{1}+1});
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
    Content = get(handles.popupDisplayOption,'String');
    SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
    if sum(strcmp(SelectedCube,{'physicalDose','effect','RBExDose','alpha','beta','RBE'})) > 0
         Suffix = '';
    else
        Idx    = find(SelectedCube == '_');
        Suffix = SelectedCube(Idx:end);
    end
    
    mPhysDose = Result.(['physicalDose' Suffix]); 
    PlotHandles{1} = plot(handles.axesFig,vX,mPhysDose(ix),'color',cColor{1,1},'LineWidth',3); hold on; 
    PlotHandles{1,2} ='physicalDose';
    ylabel(handles.axesFig,'dose in [Gy]');
    set(handles.axesFig,'FontSize',defaultFontSize);
    
    % plot counter
    Cnt=2;
    
    if isfield(Result,['RBE' Suffix])
        
        %disbale specific plots
        %DispInfo{6,2}=0;
        %DispInfo{5,2}=0;
        %DispInfo{2,2}=0;
        
        % generate two lines for ylabel
        StringYLabel1 = '\fontsize{8}{\color{red}RBE x dose [Gy(RBE)] \color{black}dose [Gy] ';
        StringYLabel2 = '';
        for i=1:1:size(DispInfo,1)
            if DispInfo{i,2} && sum(strcmp(DispInfo{i,1},{['effect' Suffix],['alpha' Suffix],['beta' Suffix]})) > 0
                %physicalDose is already plotted and RBExD vs RBE is plotted later with plotyy
                if ~strcmp(DispInfo{i,1},['RBExDose' Suffix]) &&...
                   ~strcmp(DispInfo{i,1},['RBE' Suffix]) && ...
                   ~strcmp(DispInfo{i,1},['physicalDose' Suffix])
               
                        mCube = Result.([DispInfo{i,1}]);
                        PlotHandles{Cnt,1} = plot(handles.axesFig,vX,mCube(ix),'color',cColor{1,Cnt},'LineWidth',3);hold on; 
                        PlotHandles{Cnt,2} = DispInfo{i,1};
                        StringYLabel2 = [StringYLabel2  ' \color{'  cColor{1,Cnt} '}' DispInfo{i,1} ' ['  DispInfo{i,3} ']'];
                        Cnt = Cnt+1;
                end    
            end
        end
        StringYLabel2 = [StringYLabel2 '}'];
        % always plot RBExD against RBE
        mRBExDose = Result.(['RBExDose' Suffix]);
        vBED = mRBExDose(ix);
        mRBE = Result.(['RBE' Suffix]);
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
        if strcmp(cst{i,3},'TARGET') && tmpPrior >= cst{i,5}.Priority && tmpSize<numel(cst{i,4}{1})
           linIdxTarget = unique(cst{i,4}{1});
           tmpPrior=cst{i,5}.Priority;
           tmpSize=numel(cst{i,4}{1});
           VOI = cst{i,2};
        end
    end
    
    str = sprintf('profile plot - central axis of %d beam gantry angle %d? couch angle %d?',...
        handles.selectedBeam ,pln.propStf.gantryAngles(handles.selectedBeam),pln.propStf.couchAngles(handles.selectedBeam));
    h_title = title(handles.axesFig,str,'FontSize',defaultFontSize);
    pos = get(h_title,'Position');
    set(h_title,'Position',[pos(1)-40 pos(2) pos(3)])
    
    % plot target boundaries
    mTargetCube = zeros(ct.cubeDim);
    mTargetCube(linIdxTarget) = 1;
    vProfile = mTargetCube(ix);
    WEPL_Target_Entry = vX(find(vProfile,1,'first'));
    WEPL_Target_Exit  = vX(find(vProfile,1,'last'));
    PlotHandles{Cnt,2} =[VOI ' boundary'];
    
    if ~isempty(WEPL_Target_Entry) && ~isempty(WEPL_Target_Exit)
        hold on
        PlotHandles{Cnt,1} = ...
        plot([WEPL_Target_Entry WEPL_Target_Entry],get(handles.axesFig,'YLim'),'--','Linewidth',3,'color','k');hold on
        plot([WEPL_Target_Exit WEPL_Target_Exit],get(handles.axesFig,'YLim'),'--','Linewidth',3,'color','k');hold on
      
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

zoom(handles.figure1,'reset');
axis(handles.axesFig,'tight');

if handles.rememberCurrAxes
    axis(currAxes);
end

hold(handles.axesFig,'off');

handles.cBarChanged = false;
guidata(handles.axesFig,handles);
%guidata(gcf,handles);  
if get(handles.popupTypeOfPlot,'Value')==1 
    UpdateColormapOptions(handles);
end

Update3DView(handles);

%profile off;
%profile viewer;

function Update3DView(handles)

if isfield(handles,'axesFig3D') && isfield(handles,'fig3D') && isgraphics(handles.axesFig3D) && isgraphics(handles.fig3D)
    axesFig3D = handles.axesFig3D;
    fig3D = handles.fig3D;    
else
    return
end

if handles.State == 0
    return
elseif handles.State > 0
     AllVarNames = evalin('base','who');
    if  ismember('resultGUI',AllVarNames)
        Result = evalin('base','resultGUI');
    end
    
    if  ismember('stf',AllVarNames)
        stf = evalin('base','stf');
    else
        stf = [];
    end

    ct  = evalin('base','ct');
    cst = evalin('base','cst');
    pln = evalin('base','pln');
end

oldView = get(axesFig3D,'View');

cla(axesFig3D);
%delete(allchild(axesFig3D));

%test = allchild(axesFig3D);

plane = get(handles.popupPlane,'Value');
slice = round(get(handles.sliderSlice,'Value'));
defaultFontSize = 8;

%Check if we need to precompute the surface data
if size(cst,2) < 8
    cst = matRad_computeAllVoiSurfaces(ct,cst);
    assignin('base','cst',cst);
end

set(fig3D,'Color',0.5*[1 1 1]);
set(axesFig3D,'Color',1*[0 0 0]);

%% Plot 3D structures
hold(axesFig3D,'on');
if get(handles.radiobtnContour,'Value') && handles.State>0
    voiPatches = matRad_plotVois3D(axesFig3D,ct,cst,handles.VOIPlotFlag,colorcube);
end

%% plot the CT slice
if get(handles.radiobtnCT,'Value')
    window = handles.dispWindow{2,1}; %(2 for ct)
    ctMap = matRad_getColormap(handles.ctColorMap,handles.cMapSize);
    ctHandle = matRad_plotCtSlice3D(axesFig3D,ct,1,plane,slice,ctMap,window);
end

%% plot the dose slice
if handles.State >= 1 && exist('Result','var')
    doseMap = matRad_getColormap(handles.doseColorMap,handles.cMapSize);
    doseIx  = 3;
    % if the selected display option doesn't exist then simply display
    % the first cube of the Result struct
    if ~isfield(Result,handles.SelectedDisplayOption)
        CubeNames = fieldnames(Result);
        handles.SelectedDisplayOption = CubeNames{1,1};
    end
    
    dose = Result.(handles.SelectedDisplayOption);
    
    % dose colorwash
    if ~isempty(dose) && ~isvector(dose)
        
        if isempty(handles.dispWindow{doseIx,2})
            handles.dispWindow{doseIx,2} = [min(dose(:)) max(dose(:))];   % set min and max dose values
        end
        
        if get(handles.radiobtnDose,'Value')
            [doseHandle,~,handles.dispWindow{doseIx,1}] = matRad_plotDoseSlice3D(axesFig3D,ct,dose,plane,slice,handles.CutOffLevel,handles.doseOpacity,doseMap,handles.dispWindow{doseIx,1});
        end
        if get(handles.radiobtnIsoDoseLines,'Value')
            matRad_plotIsoDoseLines3D(axesFig3D,ct,dose,handles.IsoDose.Contours,handles.IsoDose.Levels,plane,slice,doseMap,handles.dispWindow{doseIx,1},'LineWidth',1.5);
        end
    end
end

if get(handles.radiobtnPlan,'Value')
    matRad_plotPlan3D(axesFig3D,pln,stf);
end

%hLight = light('Parent',axesFig3D);
%camlight(hLight,'left');
%lighting('gouraud');

xlabel(axesFig3D,'x [voxels]','FontSize',defaultFontSize)
ylabel(axesFig3D,'y [voxels]','FontSize',defaultFontSize)
zlabel(axesFig3D,'z [voxels]','FontSize',defaultFontSize)
title(axesFig3D,'matRad 3D view');

% set axis ratio
ratios = [1 1 1]; %[1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
ratios = ratios([2 1 3]);
set(axesFig3D,'DataAspectRatioMode','manual');
set(axesFig3D,'DataAspectRatio',ratios./max(ratios));

set(axesFig3D,'Ydir','reverse');

set(axesFig3D,'view',oldView);


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
        set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(handles.plane),...
                'SliderStep',[1/(ct.cubeDim(handles.plane)-1) 1/(ct.cubeDim(handles.plane)-1)]);
        if handles.State < 3
            set(handles.sliderSlice,'Value',round(ct.cubeDim(handles.plane)/2));
        else
            pln = evalin('base','pln');
            
            if handles.plane == 1
                set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,2)/ct.resolution.x));
            elseif handles.plane == 2
                set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,1)/ct.resolution.y));
            elseif handles.plane == 3
                set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,3)/ct.resolution.z));
            end
            
        end
    end
catch
end

handles.rememberCurrAxes = false;
UpdatePlot(handles);
handles.rememberCurrAxes = true;
guidata(hObject,handles);

% --- Executes on slider movement.
function sliderSlice_Callback(~, ~, handles)
% hObject    handle to sliderSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
UpdatePlot(handles)

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
    Figures = gcf;%findobj('type','figure');
    set(Figures, 'pointer', 'watch'); 
    drawnow;
    % disable all active objects
    InterfaceObj = findobj(Figures,'Enable','on');
    set(InterfaceObj,'Enable','off');
    % wait until the table is updated
    btnTableSave_Callback([],[],handles); %We don't need it?


    % if a critical change to the cst has been made which affects the dij matrix
    if handles.DijCalcWarning == true

        choice = questdlg('Overlap priorites of OAR constraints have been edited, a new OAR VOI was added or a critical row constraint was deleted. A new Dij calculation might be necessary.', ...
        'Title','Cancel','Calculate Dij then Optimize','Optimze directly','Optimze directly');

        switch choice
            case 'Cancel'
                set(Figures, 'pointer', 'arrow');
                set(InterfaceObj,'Enable','on');
                guidata(hObject,handles);
                return;
            case 'Calculate dij again and optimize'
                handles.DijCalcWarning = false;
                btnCalcDose_Callback(hObject, eventdata, handles)      
            case 'Optimze directly'
                handles.DijCalcWarning = false;       
        end
    end
    
    pln = evalin('base','pln');
    ct  = evalin('base','ct');
    
    % optimize
    if get(handles.radiobutton3Dconf,'Value') && strcmp(handles.Modalities{get(handles.popupRadMode,'Value')},'photons')
        % conformal plan if photons and 3d conformal
        if ~matRad_checkForConnectedBixelRows(evalin('base','stf'))
            error('disconnetced dose influence data in BEV - run dose calculation again with consistent settings');
        end
        [resultGUIcurrentRun,usedOptimizer] = matRad_fluenceOptimization(matRad_collapseDij(evalin('base','dij')),evalin('base','cst'),pln);
        resultGUIcurrentRun.w = resultGUIcurrentRun.w * ones(evalin('base','dij.totalNumOfBixels'),1);
        resultGUIcurrentRun.wUnsequenced = resultGUIcurrentRun.w;
    else
        if pln.propOpt.runDAO
        if ~matRad_checkForConnectedBixelRows(evalin('base','stf'))
            error('disconnetced dose influence data in BEV - run dose calculation again with consistent settings');
        end
        end
        
        [resultGUIcurrentRun,usedOptimizer] = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),pln);
    end
    
    %if resultGUI already exists then overwrite the "standard" fields
    AllVarNames = evalin('base','who');
    if  ismember('resultGUI',AllVarNames)
        resultGUI = evalin('base','resultGUI');
        sNames = fieldnames(resultGUIcurrentRun);
        oldNames = fieldnames(resultGUI);
        if(length(oldNames) > length(sNames))
            for j = 1:length(oldNames)
            if strfind(oldNames{j}, 'beam')
                resultGUI = rmfield(resultGUI, oldNames{j});
            end
            end
        end
        for j = 1:length(sNames)
            resultGUI.(sNames{j}) = resultGUIcurrentRun.(sNames{j});
        end
    else
        resultGUI = resultGUIcurrentRun;
    end
    assignin('base','resultGUI',resultGUI);

    % set some values
    if handles.plane == 1
        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,2)/ct.resolution.x));
    elseif handles.plane == 2
        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,1)/ct.resolution.y));
    elseif handles.plane == 3
        set(handles.sliderSlice,'Value',ceil(pln.propStf.isoCenter(1,3)/ct.resolution.z));
    end

    handles.State = 3;
    handles.SelectedDisplayOptionIdx = 1;
    if strcmp(pln.radiationMode,'carbon') || (strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD'))
        handles.SelectedDisplayOption = 'RBExDose';
    else
        handles.SelectedDisplayOption = 'physicalDose';
    end
    handles.selectedBeam = 1;
    % check IPOPT status and return message for GUI user if no DAO or
    % particles
    if ~pln.propOpt.runDAO || ~strcmp(pln.radiationMode,'photons')
        CheckOptimizerStatus(usedOptimizer,'Fluence')
    end
    
catch ME
    handles = showError(handles,'OptimizeCallback: Could not optimize!',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

% perform sequencing and DAO
try
    
    %% sequencing
    if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
    %   resultGUI = matRad_xiaLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
    %       ,get(handles.editSequencingLevel,'Value'));
    %   resultGUI = matRad_engelLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
    %       ,str2double(get(handles.editSequencingLevel,'String')));
        resultGUI = matRad_siochiLeafSequencing(resultGUI,evalin('base','stf'),evalin('base','dij')...
            ,str2double(get(handles.editSequencingLevel,'String')));
        
        assignin('base','resultGUI',resultGUI);
    end
        
catch ME
    handles = showError(handles,'OptimizeCallback: Could not perform sequencing',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

try
    %% DAO
    if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
        handles = showWarning(handles,['Observe: You are running direct aperture optimization' filesep 'This is experimental code that has not been thoroughly debugged - especially in combination with constrained optimization.']);
       [resultGUI,usedOptimizer] = matRad_directApertureOptimization(evalin('base','dij'),evalin('base','cst'),...
           resultGUI.apertureInfo,resultGUI,pln);
       assignin('base','resultGUI',resultGUI);
       % check IPOPT status and return message for GUI user
       CheckOptimizerStatus(usedOptimizer,'DAO');      
    end
    
    if strcmp(pln.radiationMode,'photons') && (pln.propOpt.runSequencing || pln.propOpt.runDAO)
        matRad_visApertureInfo(resultGUI.apertureInfo);
    end
   
catch ME
    handles = showError(handles,'OptimizeCallback: Could not perform direct aperture optimization',ME); 
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    guidata(hObject,handles);
    return;
end

% change state from busy to normal
set(Figures, 'pointer', 'arrow');
set(InterfaceObj,'Enable','on');
handles.dispWindow{3,1}  = [];   % reset dose ranges
handles.dispWindow{3,2}  = [];   % reset min max dose values
handles.rememberCurrAxes = false;
handles.IsoDose.Levels   = 0;  % ensure to use default iso dose line spacing
handles.cBarChanged      = true;
    
guidata(hObject,handles);
handles = updateIsoDoseLineCache(handles);
UpdateState(handles);
UpdatePlot(handles);
handles.rememberCurrAxes = true;   
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

 % intensity plot
if get(hObject,'Value') == 1  
    
    set(handles.sliderBeamSelection,'Enable','off')
    set(handles.sliderOffset,'Enable','off')
    set(handles.popupDisplayOption,'Enable','on')
    set(handles.btnProfileType,'Enable','off');
    set(handles.popupPlane,'Enable','on');
    set(handles.radiobtnCT,'Enable','on');
    set(handles.radiobtnContour,'Enable','on');
    set(handles.radiobtnDose,'Enable','on');
    set(handles.radiobtnIsoDoseLines,'Enable','on');
    set(handles.radiobtnIsoDoseLinesLabels,'Enable','on');
    set(handles.sliderSlice,'Enable','on');
    
% profile plot
elseif get(hObject,'Value') == 2
    
    if handles.State > 0
        if length(parseStringAsNum(get(handles.editGantryAngle,'String'),true)) > 1
            
            set(handles.sliderBeamSelection,'Enable','on');
            handles.selectedBeam = 1;
            pln = evalin('base','pln');
            set(handles.sliderBeamSelection,'Min',handles.selectedBeam,'Max',pln.propStf.numOfBeams,...
                'Value',handles.selectedBeam,...
                'SliderStep',[1/(pln.propStf.numOfBeams-1) 1/(pln.propStf.numOfBeams-1)],...
                'Enable','on');

        else
            handles.selectedBeam = 1;
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
    
    
    set(handles.popupDisplayOption,'Enable','on');
    set(handles.btnProfileType,'Enable','on');
    set(handles.popupPlane,'Enable','off');
    set(handles.radiobtnCT,'Enable','off');
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

handles.cBarChanged = true;

handles.rememberCurrAxes = false;
cla(handles.axesFig,'reset');
UpdatePlot(handles);
handles.rememberCurrAxes = true;
guidata(hObject, handles);

% --- Executes on selection change in popupDisplayOption.
function popupDisplayOption_Callback(hObject, ~, handles)
content = get(hObject,'String');
handles.SelectedDisplayOption = content{get(hObject,'Value'),1};
handles.SelectedDisplayOptionIdx = get(hObject,'Value');
%handles.dispWindow{3,1} = []; handles.dispWindow{3,2} = [];

if ~isfield(handles,'colormapLocked') || ~handles.colormapLocked
    handles.dispWindow{3,1} = []; handles.dispWindow{3,2} = [];
end

handles = updateIsoDoseLineCache(handles);
handles.cBarChanged = true;
guidata(hObject, handles);
UpdatePlot(handles);
guidata(hObject, handles);

% --- Executes on slider movement.
function sliderBeamSelection_Callback(hObject, ~, handles)
% hObject    handle to sliderBeamSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles.selectedBeam = round(get(hObject,'Value'));
set(hObject, 'Value', handles.selectedBeam);
handles.rememberCurrAxes = false;
UpdatePlot(handles);
handles.rememberCurrAxes = true;
guidata(hObject,handles);

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
    
    handles.rememberCurrAxes = false;
    UpdatePlot(handles);
    handles.rememberCurrAxes = true;

    guidata(hObject, handles);
 
end



function Flag = getCstTable (handles)

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
            
            % Obj Func / constraint
            if isempty(data{j,4}) ||~isempty(strfind(data{j,4}, 'Select'))
               FlagValidParameters=false;
            else
                 NewCst{Cnt,4}(CntObjF,1).type = data{j,4};
            end
         
            % get further parameter
            if FlagValidParameters
                
              NewCst{Cnt,4}(CntObjF,1).dose       = data{j,6};
              NewCst{Cnt,4}(CntObjF,1).penalty    = data{j,5};
              NewCst{Cnt,4}(CntObjF,1).EUD        = data{j,7};
              NewCst{Cnt,4}(CntObjF,1).volume     = data{j,8};
              NewCst{Cnt,4}(CntObjF,1).robustness = data{j,9};
             
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
                   OldCst(m,6) = NewCst(n,4);
                   OldCst(m,3) = NewCst(n,2);
                   OldCst{m,5}.Priority = NewCst{n,3};
                   break;
               end 
           end

           if ~boolChanged
               OldCst{m,6}=[];
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
data{sEnd+1,4} = 'Select obj func/constraint';
data{sEnd+1,9} = 'none';

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

Placeholder = NaN;
PlaceholderRob  = 'none';

% get table data and current index of cell
if isempty(eventdata)
    data = get(handles.uiTable,'Data');
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
         cst = evalin('base','cst');
         Idx = ~cellfun('isempty',cst(:,6));
         
         if sum(strcmp(cst(Idx,2),eventdata.NewData))==0 
             handles.DijCalcWarning = true;
         end
        
        
    else
        % if table changed after a optimization was performed
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
    ObjFunction = eventdata.NewData;
else
    ObjFunction = data{eventdata.Indices(1),4};
end


if eventdata.Indices(2) == 4
    if  sum(strcmp(eventdata.NewData,{'square deviation','square underdosing','min dose constraint',...
                                      'min mean dose constraint','min DVH constraint','min DVH objective'})) > 0 
    
        if strcmp('OAR',data{eventdata.Indices(1),2})
            data{eventdata.Indices(1),4} = 'square overdosing';
            ObjFunction                  = 'square overdosing';
            if isnan(data{eventdata.Indices(1),5})
                data{eventdata.Indices(1),5} = 1;
            end
        end
        
        tmpDose = parseStringAsNum(data{eventdata.Indices(1),6},true);
        if numel(tmpDose) == 2
            data{eventdata.Indices(1),6} = num2str(tmpDose(1));
        end
    
    elseif strcmp(eventdata.NewData,'mean')
        
        if strcmp('TARGET',data{eventdata.Indices(1),2})
            data{eventdata.Indices(1),4} = 'square deviation';
        end
        
    end
end

%% set fields to NaN according to objective function

if sum(strcmp(ObjFunction, {'square underdosing','square overdosing','square deviation'})) > 0
    
    for k = [5 6]
        if isnan(data{eventdata.Indices(1),k})
             data{eventdata.Indices(1),k} = 1;
        end
    end 
    data{eventdata.Indices(1),7} = Placeholder;
    data{eventdata.Indices(1),8} = Placeholder;
    data{eventdata.Indices(1),9} = PlaceholderRob;
   
elseif strcmp(ObjFunction,'mean')
    
        if isnan(data{eventdata.Indices(1),5})
                 data{eventdata.Indices(1),5} = 1;
        end
        data{eventdata.Indices(1),6} = Placeholder;    
        data{eventdata.Indices(1),7} = Placeholder;   
        data{eventdata.Indices(1),8} = Placeholder;
        data{eventdata.Indices(1),9} = PlaceholderRob;

elseif strcmp(ObjFunction,'EUD')
    
        for k = [5 7]
            if isnan(data{eventdata.Indices(1),k})
                 data{eventdata.Indices(1),k} = 1;
            end
        end 
       data{eventdata.Indices(1),6} = Placeholder; 
       data{eventdata.Indices(1),8} = Placeholder;
       data{eventdata.Indices(1),9} = PlaceholderRob;
       
elseif sum(strcmp(ObjFunction,{'min dose constraint','max dose constraint'...
                               'min mean dose constraint','max mean dose constraint'}))> 0
         
         if isnan(data{eventdata.Indices(1),6})
                 data{eventdata.Indices(1),6} = 1;
         end
         data{eventdata.Indices(1),5} = Placeholder;
         data{eventdata.Indices(1),7} = Placeholder;
         data{eventdata.Indices(1),8} = Placeholder;
         data{eventdata.Indices(1),9} = PlaceholderRob;
         
elseif sum(strcmp(ObjFunction,{'min EUD constraint','max EUD constraint'}) ) > 0
        
        if isnan(data{eventdata.Indices(1),7})
             data{eventdata.Indices(1),7} = 1;
        end
        data{eventdata.Indices(1),5} = Placeholder;
        data{eventdata.Indices(1),6} = Placeholder;
        data{eventdata.Indices(1),8} = Placeholder;
        data{eventdata.Indices(1),9} = PlaceholderRob;
        
elseif sum(strcmp(ObjFunction,{'min DVH constraint','max DVH constraint'}) ) > 0
        
        for k = [6 8]
            if isnan(data{eventdata.Indices(1),k})
                 data{eventdata.Indices(1),k} = 1;
            end
        end    
        data{eventdata.Indices(1),5} = Placeholder;
        data{eventdata.Indices(1),7} = Placeholder;
        data{eventdata.Indices(1),9} = PlaceholderRob;

elseif sum(strcmp(ObjFunction,{'min DVH objective','max DVH objective'}) ) > 0
        
    for k = [5 6 8]
        if isnan(data{eventdata.Indices(1),k})
             data{eventdata.Indices(1),k} = 1;
        end
    end
    data{eventdata.Indices(1),7} = Placeholder;
    data{eventdata.Indices(1),9} = PlaceholderRob;
    
end
    
%% check if input is a valid
%check if overlap, penalty and and parameters are numbers
if (eventdata.Indices(2) == 3  || eventdata.Indices(2) == 5 || eventdata.Indices(2) == 6 || eventdata.Indices(2) == 7 || eventdata.Indices(2) == 8) ...
        && ~isempty(eventdata.NewData)
    if CheckValidity(eventdata.NewData) == false
            data{eventdata.Indices(1),eventdata.Indices(2)} = eventdata.PreviousData;
    end
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

% enables/ disables buttons according to the current state      
function UpdateState(handles)

if handles.State > 0
    pln = evalin('base','pln');

    if strcmp(pln.radiationMode,'carbon')
        set(handles.popMenuBioOpt,'Enable','on');
        set(handles.btnSetTissue,'Enable','on');
    elseif strcmp(pln.radiationMode,'protons')
        set(handles.popMenuBioOpt,'Enable','on');
        set(handles.btnSetTissue,'Enable','off');
    else
        set(handles.popMenuBioOpt,'Enable','off');
        set(handles.btnSetTissue,'Enable','off'); 
    end
    
    cMapControls = allchild(handles.uipanel_colormapOptions);
    for runHandles = cMapControls
        set(runHandles,'Enable','on');
    end
end 

if handles.cubeHUavailable
    cMapOptionsSelectList = {'None','CT (HU)','Result (i.e. dose)'};
    set(handles.popupmenu_windowPreset,'Visible','on');
    set(handles.text_windowPreset,'String','Window Preset');
else
    cMapOptionsSelectList = {'None','CT (ED)','Result (i.e. dose)'};
    set(handles.popupmenu_windowPreset,'Visible','off');
    set(handles.text_windowPreset,'String','No available Window Presets');
end
handles.cBarChanged = true;

 switch handles.State
     
     case 0
      
      set(handles.txtInfo,'String','no data loaded');
      set(handles.btnCalcDose,'Enable','off');
      set(handles.btnOptimize ,'Enable','off');
      set(handles.pushbutton_recalc,'Enable','off');
      set(handles.btnSaveToGUI,'Enable','off');
      set(handles.btnDVH,'Enable','off'); 
      set(handles.importDoseButton,'Enable','off');
      set(handles.btn_export,'Enable','off');
      set(handles.btn3Dview,'Enable','off');
      
      cMapControls = allchild(handles.uipanel_colormapOptions);
      for runHandles = cMapControls
          set(runHandles,'Enable','off');
      end
      
      set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList{1})
      set(handles.popupmenu_chooseColorData,'Value',1);
      
     case 1
     
      set(handles.txtInfo,'String','ready for dose calculation');
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','off');
      set(handles.pushbutton_recalc,'Enable','off');
      set(handles.btnSaveToGUI,'Enable','off');
      set(handles.btnDVH,'Enable','off');
      set(handles.importDoseButton,'Enable','off');
      set(handles.btn_export,'Enable','on');
      set(handles.btn3Dview,'Enable','on');
      
      set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList(1:2))
      set(handles.popupmenu_chooseColorData,'Value',2);
      AllVarNames = evalin('base','who');
      if ~isempty(AllVarNames)
            if  ismember('resultGUI',AllVarNames)
              set(handles.pushbutton_recalc,'Enable','on');
              set(handles.btnSaveToGUI,'Enable','on');
              set(handles.btnDVH,'Enable','on');
              set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList(1:3))
              set(handles.popupmenu_chooseColorData,'Value',3);
            end
      end

     case 2
    
      set(handles.txtInfo,'String','ready for optimization');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.pushbutton_recalc,'Enable','off');
      set(handles.btnSaveToGUI,'Enable','off');
      set(handles.btnDVH,'Enable','off');
      set(handles.importDoseButton,'Enable','off');
      set(handles.btn_export,'Enable','on');
      set(handles.btn3Dview,'Enable','on');
      set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList(1:2))
      set(handles.popupmenu_chooseColorData,'Value',2);
      AllVarNames = evalin('base','who');
      
      if ~isempty(AllVarNames)
            if  ismember('resultGUI',AllVarNames)
              set(handles.pushbutton_recalc,'Enable','on');
              set(handles.btnSaveToGUI,'Enable','on');
              set(handles.btnDVH,'Enable','on');
              set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList(1:3))
              set(handles.popupmenu_chooseColorData,'Value',3);
            end
      end
      
     case 3
      set(handles.txtInfo,'String','plan is optimized');   
      set(handles.btnCalcDose,'Enable','on');
      set(handles.btnOptimize ,'Enable','on');
      set(handles.pushbutton_recalc,'Enable','on');
      set(handles.btnSaveToGUI,'Enable','on');
      set(handles.btnDVH,'Enable','on');
      set(handles.btn_export,'Enable','on');
      set(handles.btn3Dview,'Enable','on');
      % resultGUI struct needs to be available to import dose
      % otherwise inconsistent states can be achieved
      set(handles.importDoseButton,'Enable','on');
      set(handles.popupmenu_chooseColorData,'String',cMapOptionsSelectList(1:3))
      set(handles.popupmenu_chooseColorData,'Value',3);
 end

guidata(handles.figure1,handles); 
 
% fill GUI elements with plan information
function setPln(handles)
pln = evalin('base','pln');
% sanity check of isoCenter
if size(pln.propStf.isoCenter,1) ~= pln.propStf.numOfBeams && size(pln.propStf.isoCenter,1) == 1
  pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * pln.propStf.isoCenter(1,:);
elseif size(pln.propStf.isoCenter,1) ~= pln.propStf.numOfBeams && size(pln.propStf.isoCenter,1) ~= 1
  error('Isocenter in plan file are incosistent.');
end
set(handles.editBixelWidth,'String',num2str(pln.propStf.bixelWidth));
set(handles.editFraction,'String',num2str(pln.numOfFractions));

if isfield(pln.propStf,'isoCenter')
    if size(unique(pln.propStf.isoCenter,'rows'),1) == 1
        set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.propStf.isoCenter(1,:)*10))./10), '\s+', ' '));
        set(handles.editIsoCenter,'Enable','on');
        set(handles.checkIsoCenter,'Enable','on');
    else
        set(handles.editIsoCenter,'String','multiple isoCenter');
        set(handles.editIsoCenter,'Enable','off');
        set(handles.checkIsoCenter,'Value',0);
        set(handles.checkIsoCenter,'Enable','off');
    end
end
set(handles.editGantryAngle,'String',num2str((pln.propStf.gantryAngles)));
set(handles.editCouchAngle,'String',num2str((pln.propStf.couchAngles)));
set(handles.popupRadMode,'Value',find(strcmp(get(handles.popupRadMode,'String'),pln.radiationMode)));
set(handles.popUpMachine,'Value',find(strcmp(get(handles.popUpMachine,'String'),pln.machine)));

if ~strcmp(pln.propOpt.bioOptimization,'none')  
    set(handles.popMenuBioOpt,'Enable','on');
    contentPopUp = get(handles.popMenuBioOpt,'String');
    ix = find(strcmp(pln.propOpt.bioOptimization,contentPopUp));
    set(handles.popMenuBioOpt,'Value',ix);
    set(handles.btnSetTissue,'Enable','on');
else
    set(handles.popMenuBioOpt,'Enable','off');
    set(handles.btnSetTissue,'Enable','off');
end
%% enable sequencing and DAO button if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runSequencing
    set(handles.btnRunSequencing,'Enable','on');
    set(handles.btnRunSequencing,'Value',1);
elseif strcmp(pln.radiationMode,'photons') && ~pln.propOpt.runSequencing
    set(handles.btnRunSequencing,'Enable','on');
    set(handles.btnRunSequencing,'Value',0);
else
    set(handles.btnRunSequencing,'Enable','off');
end
%% enable DAO button if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons') && pln.propOpt.runDAO
    set(handles.btnRunDAO,'Enable','on');
    set(handles.btnRunDAO,'Value',1);
elseif strcmp(pln.radiationMode,'photons') && ~pln.propOpt.runDAO
    set(handles.btnRunDAO,'Enable','on');
    set(handles.btnRunDAO,'Value',0);
else
    set(handles.btnRunDAO,'Enable','off');
end
%% enable stratification level input if radiation mode is set to photons
if strcmp(pln.radiationMode,'photons')
    set(handles.txtSequencing,'Enable','on');
    set(handles.radiobutton3Dconf,'Enable','on');
    set(handles.editSequencingLevel,'Enable','on');
else
    set(handles.txtSequencing,'Enable','off');
    set(handles.radiobutton3Dconf,'Enable','off');
    set(handles.editSequencingLevel,'Enable','off');
end

% --- Executes on button press in btnTableSave.
function btnTableSave_Callback(~, ~, handles)
% hObject    handle to btnTableSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%getCstTable(handles);
if get(handles.checkIsoCenter,'Value')
    pln = evalin('base','pln'); 
    pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct')); 
    set(handles.editIsoCenter,'String',regexprep(num2str((round(pln.propStf.isoCenter(1,:) * 10))./10), '\s+', ' '));
    assignin('base','pln',pln);
end
getPlnFromGUI(handles);

% --- Executes on selection change in listBoxCmd.
function listBoxCmd_Callback(hObject, ~, ~)
numLines = size(get(hObject,'String'),1);
set(hObject, 'ListboxTop', numLines);

% --- Executes on slider movement.
function sliderOffset_Callback(hObject, ~, handles)
handles.profileOffset = get(hObject,'Value');
UpdatePlot(handles);

 
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check validity of input for cst
function flagValidity = CheckValidity(Val) 
      
flagValidity = true;

if ischar(Val)
    Val = str2num(Val);
end 

if length(Val) > 2
    warndlg('invalid input!');
end

if isempty(Val)
   warndlg('Input not a number!');
   flagValidity = false;        
end

if any(Val < 0) 
   warndlg('Input not a positive number!');
   flagValidity = false;  
end

% return IPOPT status as message box
function CheckOptimizerStatus(usedOptimizer,OptCase) 
      
[statusmsg,statusflag] = usedOptimizer.GetStatus();
    
if statusflag == 0 || statusflag == 1
    status = 'none';
else
    status = 'warn';
end

msgbox(['Optimizer finished with status ' num2str(statusflag) ' (' statusmsg ')'],'Optimizer',status,'modal');

% get pln file form GUI     
function getPlnFromGUI(handles)

% evalin pln (if existant) in order to decide whether isoCenter should be calculated
% automatically
if evalin('base','exist(''pln'',''var'')')
    pln = evalin('base','pln');
end

pln.propStf.bixelWidth      = parseStringAsNum(get(handles.editBixelWidth,'String'),false); % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles    = parseStringAsNum(get(handles.editGantryAngle,'String'),true); % [???]
pln.propStf.couchAngles     = parseStringAsNum(get(handles.editCouchAngle,'String'),true); % [???]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
try
    ct = evalin('base','ct');
    pln.numOfVoxels     = prod(ct.cubeDim);
    pln.voxelDimensions = ct.cubeDim;
catch
end
pln.numOfFractions  = parseStringAsNum(get(handles.editFraction,'String'),false);
contents            = get(handles.popupRadMode,'String'); 
pln.radiationMode   = contents{get(handles.popupRadMode,'Value')}; % either photons / protons / carbon
contents            = get(handles.popUpMachine,'String'); 
pln.machine         = contents{get(handles.popUpMachine,'Value')}; 

if (~strcmp(pln.radiationMode,'photons'))
    contentBioOpt = get(handles.popMenuBioOpt,'String');
    pln.propOpt.bioOptimization = contentBioOpt{get(handles.popMenuBioOpt,'Value'),:};
else
    pln.propOpt.bioOptimization = 'none';
end

pln.propOpt.runSequencing = logical(get(handles.btnRunSequencing,'Value'));
pln.propOpt.runDAO = logical(get(handles.btnRunDAO,'Value'));

try
    cst = evalin('base','cst');
    if (sum(strcmp('TARGET',cst(:,3))) > 0 && get(handles.checkIsoCenter,'Value')) || ...
            (sum(strcmp('TARGET',cst(:,3))) > 0 && ~isfield(pln.propStf,'isoCenter'))
       pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct);
       set(handles.checkIsoCenter,'Value',1);
    else
        if ~strcmp(get(handles.editIsoCenter,'String'),'multiple isoCenter')
            pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1) * str2num(get(handles.editIsoCenter,'String'));
        end
    end
catch
    warning('couldnt set isocenter in getPln function')
end

handles.pln = pln;
assignin('base','pln',pln);

% parsing a string as number array
function number = parseStringAsNum(stringIn,isVector)
if isnumeric(stringIn)
    number = stringIn;
else
    number = str2num(stringIn);
    if isempty(number) || length(number) > 1 && ~isVector
        warndlg(['could not parse all parameters (pln, optimization parameter)']); 
        number = NaN;
    elseif isVector && iscolumn(number)
        number = number';
    end
end


% show error
function handles = showError(handles,Message,ME)

if nargin == 3
    %Add exception message
    if isfield(handles,'devMode') && handles.devMode
        meType = 'extended';
    else 
        meType = 'basic';
    end
    Message = {Message,ME.getReport(meType,'hyperlinks','off')};    
end

if isfield(handles,'ErrorDlg')
    if ishandle(handles.ErrorDlg)
        close(handles.ErrorDlg);
    end
end
handles.ErrorDlg = errordlg(Message);

% show warning
function handles = showWarning(handles,Message,ME)

if nargin == 3
    %Add exception message
    if isfield(handles,'devMode') && handles.devMode
        meType = 'extended';
    else 
        meType = 'basic';
    end
    Message = {Message,ME.getReport(meType,'hyperlinks','off')};    
end

if isfield(handles,'WarnDlg')
    if ishandle(handles.WarnDlg)
        close(handles.WarnDlg);
    end
end
handles.WarnDlg = warndlg(Message);

% check for valid machine data input file
function flag = checkRadiationComposition(handles)
flag = true;
contents = cellstr(get(handles.popUpMachine,'String'));
Machine = contents{get(handles.popUpMachine,'Value')};
contents = cellstr(get(handles.popupRadMode,'String'));
radMod = contents{get(handles.popupRadMode,'Value')};

if isdeployed
    FoundFile = dir([ctfroot filesep 'matRad' filesep radMod '_' Machine '.mat']);
else
    FoundFile = dir([fileparts(mfilename('fullpath')) filesep radMod '_' Machine '.mat']);    
end
if isempty(FoundFile)
    warndlg(['No base data available for machine: ' Machine]);
    flag = false;
end

function matRadScrollWheelFcn(src,event)

% get handles
handles = guidata(src);

% compute new slice
currSlice = round(get(handles.sliderSlice,'Value'));
newSlice  = currSlice - event.VerticalScrollCount;

% project to allowed set
newSlice = min(newSlice,get(handles.sliderSlice,'Max'));
newSlice = max(newSlice,get(handles.sliderSlice,'Min'));

% update slider
set(handles.sliderSlice,'Value',newSlice);

% update handles object
guidata(src,handles);

% update plot
UpdatePlot(handles);




%% CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% button: show DVH
function btnDVH_Callback(~, ~, handles)

resultGUI = evalin('base','resultGUI');
Content = get(handles.popupDisplayOption,'String');
SelectedCube = Content{get(handles.popupDisplayOption,'Value')};

pln = evalin('base','pln');
resultGUI_SelectedCube.physicalDose = resultGUI.(SelectedCube);

if ~strcmp(pln.propOpt.bioOptimization,'none')

    %check if one of the default fields is selected
    if sum(strcmp(SelectedCube,{'physicalDose','effect','RBE,','RBExDose','alpha','beta'})) > 0
        resultGUI_SelectedCube.physicalDose = resultGUI.physicalDose;
        resultGUI_SelectedCube.RBExDose     = resultGUI.RBExDose;
    else
        Idx    = find(SelectedCube == '_');
        SelectedSuffix = SelectedCube(Idx(1):end);
        resultGUI_SelectedCube.physicalDose = resultGUI.(['physicalDose' SelectedSuffix]);
        resultGUI_SelectedCube.RBExDose     = resultGUI.(['RBExDose' SelectedSuffix]);
    end
end

%adapt visibilty
cst = evalin('base','cst');
for i = 1:size(cst,1)
    cst{i,5}.Visible = handles.VOIPlotFlag(i);
end

matRad_indicatorWrapper(cst,pln,resultGUI_SelectedCube);

assignin('base','cst',cst);

% radio button: plot isolines labels
function radiobtnIsoDoseLinesLabels_Callback(~, ~, handles)
UpdatePlot(handles);

% button: refresh
function btnRefresh_Callback(hObject, ~, handles)

handles = resetGUI(hObject, handles);

%% parse variables from base workspace
AllVarNames = evalin('base','who');
handles.AllVarNames = AllVarNames;
try
    if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
        ct  = evalin('base','ct');
        cst = evalin('base','cst');
        %cst = setCstTable(handles,cst);
        generateCstTable(handles,cst);
        handles.State = 1;
        cst = matRad_computeVoiContoursWrapper(cst,ct);
        assignin('base','cst',cst);
    elseif ismember('ct',AllVarNames) &&  ~ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find cst file');
    elseif ~ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
         handles = showError(handles,'GUI OpeningFunc: could not find ct file');
    end
catch  
   handles = showError(handles,'GUI OpeningFunc: Could not load ct and cst file');
end

if ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
    handles = reloadGUI(hObject, handles, ct, cst);
else
    handles = reloadGUI(hObject, handles);
end
guidata(hObject, handles);


% text box: # fractions
function editFraction_Callback(hObject, ~, handles)
getPlnFromGUI(handles);
guidata(hObject,handles);

% text box: stratification levels
function editSequencingLevel_Callback(~, ~, ~)

% text box: isoCenter in [mm]
function editIsoCenter_Callback(hObject, ~, handles)

pln = evalin('base','pln');
tmpIsoCenter = str2num(get(hObject,'String'));

if length(tmpIsoCenter) == 3
    if sum(any(unique(pln.propStf.isoCenter,'rows')~=tmpIsoCenter))
        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
        handles.State = 1;
        UpdateState(handles);
    end
else
    handles = showError(handles,'EditIsoCenterCallback: Could not set iso center');
end

assignin('base','pln',pln);
guidata(hObject,handles);

% check box: iso center auto
function checkIsoCenter_Callback(hObject, ~, handles)

W = evalin('base','whos');
doesPlnExist = ismember('pln',{W(:).name});

if get(hObject,'Value') && doesPlnExist
    pln = evalin('base','pln');
    if ~isfield(pln.propStf,'isoCenter')
        pln.propStf.isoCenter = NaN;
    end
    tmpIsoCenter = matRad_getIsoCenter(evalin('base','cst'),evalin('base','ct'));
    if ~isequal(tmpIsoCenter,pln.propStf.isoCenter)
        pln.propStf.isoCenter = ones(pln.propStf.numOfBeams,1)*tmpIsoCenter;
        handles.State = 1;
        UpdateState(handles);
    end
    set(handles.editIsoCenter,'String',regexprep(num2str((round(tmpIsoCenter*10))./10), '\s+', ' '));
    set(handles.editIsoCenter,'Enable','off')
    assignin('base','pln',pln);
else
    set(handles.editIsoCenter,'Enable','on')
end

% radio button: run sequencing
function btnRunSequencing_Callback(~, ~, handles)
getPlnFromGUI(handles);

% radio button: run direct aperture optimization
function btnRunDAO_Callback(~, ~, handles)
getPlnFromGUI(handles);

% button: set iso dose levels
function btnSetIsoDoseLevels_Callback(hObject, ~, handles)
prompt = {['Enter iso dose levels in [Gy]. Enter space-separated numbers, e.g. 1.5 2 3 4.98. Enter 0 to use default values']};
if isequal(handles.IsoDose.Levels,0) || ~isvector(handles.IsoDose.Levels) || any(~isnumeric(handles.IsoDose.Levels)) || any(isnan(handles.IsoDose.Levels))
    defaultLine = {'1 2 3 '};
else
    if isrow(handles.IsoDose.Levels)
        defaultLine = cellstr(num2str(handles.IsoDose.Levels,'%.2g '));
    else 
        defaultLine = cellstr(num2str(handles.IsoDose.Levels','%.2g '));
    end
end

try
    Input = inputdlg(prompt,'Set iso dose levels ', [1 70],defaultLine);
    if ~isempty(Input)
        handles.IsoDose.Levels = (sort(str2num(Input{1})));
        if length(handles.IsoDose.Levels) == 1 && (handles.IsoDose.Levels(1) ~= 0)
            handles.IsoDose.Levels = [handles.IsoDose.Levels handles.IsoDose.Levels];
        end
        handles.IsoDose.NewIsoDoseFlag = true;
    end
catch
    warning('Couldnt parse iso dose levels - using default values');
    handles.IsoDose.Levels = 0;
end
handles = updateIsoDoseLineCache(handles);
handles.IsoDose.NewIsoDoseFlag = false;
UpdatePlot(handles);
guidata(hObject,handles);


% popup menu: machine
function popUpMachine_Callback(hObject, ~, handles)
contents = cellstr(get(hObject,'String'));
checkRadiationComposition(handles);
if handles.State > 0
    pln = evalin('base','pln');
    if handles.State > 0 && ~strcmp(contents(get(hObject,'Value')),pln.machine)
        handles.State = 1;
        UpdateState(handles);
        guidata(hObject,handles);
    end
   getPlnFromGUI(handles);
end

% toolbar load button
function toolbarLoad_ClickedCallback(hObject, eventdata, handles)
btnLoadMat_Callback(hObject, eventdata, handles);

% toolbar save button
function toolbarSave_ClickedCallback(hObject, eventdata, handles)

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

% button: about
function btnAbout_Callback(hObject, eventdata, handles)

msgbox({'https://github.com/e0404/matRad/' 'email: matrad@dkfz.de'},'About');

% button: close
function figure1_CloseRequestFcn(hObject, ~, ~)
set(0,'DefaultUicontrolBackgroundColor',[0.5 0.5 0.5]);     
selection = questdlg('Do you really want to close matRad?',...
                     'Close matRad',...
                     'Yes','No','Yes');

%BackgroundColor',[0.5 0.5 0.5]
 switch selection
   case 'Yes',
     delete(hObject);
   case 'No'
     return
 end

% --- Executes on button press in pushbutton_recalc.
function pushbutton_recalc_Callback(hObject, ~, handles)

% recalculation only makes sense if ...
if evalin('base','exist(''pln'',''var'')') && ...
   evalin('base','exist(''stf'',''var'')') && ...
   evalin('base','exist(''ct'',''var'')') && ...
   evalin('base','exist(''cst'',''var'')') && ...
   evalin('base','exist(''resultGUI'',''var'')')

try

    % indicate that matRad is busy
    % change mouse pointer to hour glass 
    Figures = gcf;%findobj('type','figure');
    set(Figures, 'pointer', 'watch'); 
    drawnow;
    % disable all active objects
    InterfaceObj = findobj(Figures,'Enable','on');
    set(InterfaceObj,'Enable','off');
    
    % get all data from workspace
    pln       = evalin('base','pln');
    stf       = evalin('base','stf');
    ct        = evalin('base','ct');
    cst       = evalin('base','cst');
    resultGUI = evalin('base','resultGUI');
    
    % get weights of the selected cube
    Content = get(handles.popupDisplayOption,'String');
    SelectedCube = Content{get(handles.popupDisplayOption,'Value')};
    Suffix = strsplit(SelectedCube,'_');
    if length(Suffix)>1
        Suffix = ['_' Suffix{2}];
    else
        Suffix = '';
    end
    
    if sum([stf.totalNumOfBixels]) ~= length(resultGUI.(['w' Suffix]))
        warndlg('weight vector does not corresponding to current steering file');
        return
    end
    
    % change isocenter if that was changed and do _not_ recreate steering
    % information
    for i = 1:numel(pln.propStf.gantryAngles)
        stf(i).isoCenter = pln.propStf.isoCenter(i,:);
    end

    % recalculate influence matrix
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(ct,stf,pln,cst);
    end

    % recalculate cubes in resultGUI
    resultGUIreCalc = matRad_calcCubes(resultGUI.(['w' Suffix]),dij,cst);
    
    % delete old variables to avoid confusion
    if isfield(resultGUI,'effect')
        resultGUI = rmfield(resultGUI,'effect');
        resultGUI = rmfield(resultGUI,'RBExDose'); 
        resultGUI = rmfield(resultGUI,'RBE'); 
        resultGUI = rmfield(resultGUI,'alpha'); 
        resultGUI = rmfield(resultGUI,'beta');
    end
    
    % overwrite the "standard" fields
    sNames = fieldnames(resultGUIreCalc);
    for j = 1:length(sNames)
        resultGUI.(sNames{j}) = resultGUIreCalc.(sNames{j});
    end
    
    % assign results to base worksapce
    assignin('base','dij',dij);
    assignin('base','resultGUI',resultGUI);
    
    handles.State = 3;
    
    % show physicalDose of newly computed state
    handles.SelectedDisplayOption = 'physicalDose';
    set(handles.popupDisplayOption,'Value',find(strcmp('physicalDose',Content)));
    
    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');
    
    handles.cBarChanged = true;
    
    handles = updateIsoDoseLineCache(handles);

    UpdateState(handles);
    
    handles.rememberCurrAxes = false;
    UpdatePlot(handles);
    handles.rememberCurrAxes = true;   

    guidata(hObject,handles);

catch ME
    handles = showError(handles,'CalcDoseCallback: Error in dose recalculation!',ME); 

    % change state from busy to normal
    set(Figures, 'pointer', 'arrow');
    set(InterfaceObj,'Enable','on');

    guidata(hObject,handles);
    return;

end
    
end


% --- Executes on button press in btnSetTissue.
function btnSetTissue_Callback(hObject, ~, handles)

%check if patient is loaded
if handles.State == 0
    return
end

%parse variables from base-workspace
cst = evalin('base','cst');
pln = evalin('base','pln');

fileName = [pln.radiationMode '_' pln.machine];
load(fileName);

% check for available cell types characterized by alphaX and betaX 
for i = 1:size(machine.data(1).alphaX,2)
    CellType{i} = [num2str(machine.data(1).alphaX(i)) ' ' num2str(machine.data(1).betaX(i))];
end

%fill table data array
for i = 1:size(cst,1)
    data{i,1} = cst{i,2};
    data{i,2} = [num2str(cst{i,5}.alphaX) ' ' num2str(cst{i,5}.betaX)];
    data{i,3} = (cst{i,5}.alphaX / cst{i,5}.betaX );
end

Width  = 400;
Height = 200 + 20*size(data,1);
ScreenSize = get(0,'ScreenSize');
% show "set tissue parameter" window
figHandles = get(0,'Children');
if ~isempty(figHandles)
    IdxHandle = strcmp(get(figHandles,'Name'),'Set Tissue Parameters');
else
    IdxHandle = [];
end

%check if window is already exists 
if any(IdxHandle)
    IdxTable = find(strcmp({figHandles(IdxHandle).Children.Type},'uitable'));
    set(figHandles(IdxHandle).Children(IdxTable), 'Data', []);
    figTissue = figHandles(IdxHandle);
    %set focus
    figure(figTissue);
else
    figTissue = figure('Name','Set Tissue Parameters','Color',[.5 .5 .5],'NumberTitle','off','Position',...
        [ceil(ScreenSize(3)/2) ceil(ScreenSize(4)/2) Width Height]);
end

% define the tissue parameter table
cNames = {'VOI','alphaX betaX','alpha beta ratio'};
columnformat = {'char',CellType,'numeric'};

tissueTable = uitable('Parent', figTissue,'Data', data,'ColumnEditable',[false true false],...
                      'ColumnName',cNames, 'ColumnFormat',columnformat,'Position',[50 150 10 10]);
set(tissueTable,'CellEditCallback',@tissueTable_CellEditCallback);
% set width and height
currTablePos = get(tissueTable,'Position');
currTableExt = get(tissueTable,'Extent');
currTablePos(3) = currTableExt(3);
currTablePos(4) = currTableExt(4);
set(tissueTable,'Position',currTablePos);

% define two buttons with callbacks
uicontrol('Parent', figTissue,'Style', 'pushbutton', 'String', 'Save&Close',...
        'Position', [Width-(0.25*Width) 0.1 * Height 70 30],...
        'Callback', @(hpb,eventdata)SaveTissueParameters(hpb,eventdata,handles));
    
uicontrol('Parent', figTissue,'Style', 'pushbutton', 'String', 'Cancel&Close',...
        'Position', [Width-(0.5*Width) 0.1 * Height 80 30],...
        'Callback', 'close');    
    
guidata(hObject,handles); 
UpdateState(handles);
    
    
function SaveTissueParameters(~, ~, handles) 
cst = evalin('base','cst');    
% get handle to uiTable
figHandles = get(0,'Children');
IdxHandle  = find(strcmp(get(figHandles,'Name'),'Set Tissue Parameters'));
% find table in window

figHandleChildren = get(figHandles(IdxHandle),'Children');
IdxTable   = find(strcmp(get(figHandleChildren,'Type'),'uitable'));
uiTable    = figHandleChildren(IdxTable);
% retrieve data from uitable
data       = get(uiTable,'data');

for i = 1:size(cst,1)
   for j = 1:size(data,1)
       if strcmp(cst{i,2},data{j,1})
         alphaXBetaX = str2num(data{j,2});
         cst{i,5}.alphaX = alphaXBetaX(1);
         cst{i,5}.betaX  = alphaXBetaX(2);
       end
   end
end
assignin('base','cst',cst);
close
handles.State = 2;
UpdateState(handles);
 

        
function tissueTable_CellEditCallback(hObject, eventdata, ~) 
if eventdata.Indices(2) == 2
   alphaXBetaX = str2num(eventdata.NewData);
   data = get(hObject,'Data');
   data{eventdata.Indices(1),3} = alphaXBetaX(1)/alphaXBetaX(2);
   set(hObject,'Data',data);
end

% --- Executes on button press in btnSaveToGUI.
function btnSaveToGUI_Callback(hObject, ~, handles)

Width  = 400;
Height = 200;
ScreenSize = get(0,'ScreenSize');

% show "Provide result name" window
figHandles = get(0,'Children');
if ~isempty(figHandles)
    IdxHandle = strcmp(get(figHandles,'Name'),'Provide result name');
else
    IdxHandle = [];
end

%check if window is already exists 
if any(IdxHandle)
    figDialog = figHandles(IdxHandle);
    %set focus
    figure(figDialog);
else
    figDialog = dialog('Position',[ceil(ScreenSize(3)/2) ceil(ScreenSize(4)/2) Width Height],'Name','Provide result name','Color',[0.5 0.5 0.5]);
    
    uicontrol('Parent',figDialog,...
              'Style','text',...
              'Position',[20 Height - (0.35*Height) 350 60],...
              'String','Please provide a decriptive name for your optimization result:','FontSize',10,'BackgroundColor',[0.5 0.5 0.5]);
    
    uicontrol('Parent',figDialog,...
              'Style','edit',...
              'Position',[30 60 350 60],... 
              'String','Please enter name here...','FontSize',10,'BackgroundColor',[0.55 0.55 0.55]);
         
    uicontrol('Parent', figDialog,'Style', 'pushbutton', 'String', 'Save','FontSize',10,...
              'Position', [0.42*Width 0.1 * Height 70 30],...
              'Callback', @(hpb,eventdata)SaveResultToGUI(hpb,eventdata,guidata(hpb)));        
end

uiwait(figDialog);
guidata(hObject, handles);
UpdateState(handles)
UpdatePlot(handles)


function SaveResultToGUI(~, ~, ~)
AllFigHandles = get(0,'Children');    
ixHandle      = strcmp(get(AllFigHandles,'Name'),'Provide result name');
uiEdit        = get(AllFigHandles(ixHandle),'Children');

if strcmp(get(uiEdit(2),'String'),'Please enter name here...')
  
    formatOut = 'mmddyyHHMM';
    Suffix = ['_' datestr(now,formatOut)];
else
    % delete special characters
    Suffix = get(uiEdit(2),'String');
    logIx = isstrprop(Suffix,'alphanum');
    Suffix = ['_' Suffix(logIx)];
end

pln       = evalin('base','pln');
resultGUI = evalin('base','resultGUI');

if isfield(resultGUI,'physicalDose')
    resultGUI.(['physicalDose' Suffix])  = resultGUI.physicalDose; 
end
if isfield(resultGUI,'w')
    resultGUI.(['w' Suffix])             = resultGUI.w;
end


if ~strcmp(pln.propOpt.bioOptimization,'none')
    
    if isfield(resultGUI,'RBExDose')
         resultGUI.(['RBExDose' Suffix]) = resultGUI.RBExDose; 
    end
    
    if strcmp(pln.radiationMode,'carbon') == 1 
        if isfield(resultGUI,'effect')
            resultGUI.(['effect' Suffix])= resultGUI.effect; 
        end

        if isfield(resultGUI,'RBE')
            resultGUI.(['RBE' Suffix]) = resultGUI.RBE;
        end
        if isfield(resultGUI,'alpha')
            resultGUI.(['alpha' Suffix]) = resultGUI.alpha;
        end
        if isfield(resultGUI,'beta')
            resultGUI.(['beta' Suffix]) = resultGUI.beta;
        end
    end
end

close(AllFigHandles(ixHandle));
assignin('base','resultGUI',resultGUI);

% precompute contours of VOIs
function cst = precomputeContours(ct,cst)
mask = zeros(ct.cubeDim); % create zero cube with same dimeonsions like dose cube
for s = 1:size(cst,1)
    cst{s,7} = cell(max(ct.cubeDim(:)),3);
    mask(:) = 0;
    mask(cst{s,4}{1}) = 1;    
    for slice = 1:ct.cubeDim(1)
        if sum(sum(mask(slice,:,:))) > 0
             cst{s,7}{slice,1} = contourc(squeeze(mask(slice,:,:)),.5*[1 1]);
        end
    end
    for slice = 1:ct.cubeDim(2)
        if sum(sum(mask(:,slice,:))) > 0
             cst{s,7}{slice,2} = contourc(squeeze(mask(:,slice,:)),.5*[1 1]);
        end
    end
    for slice = 1:ct.cubeDim(3)
        if sum(sum(mask(:,:,slice))) > 0
             cst{s,7}{slice,3} = contourc(squeeze(mask(:,:,slice)),.5*[1 1]);
        end
    end
end

%Update IsodoseLines
function handles = updateIsoDoseLineCache(handles)
resultGUI = evalin('base','resultGUI');
% select first cube if selected option does not exist
if ~isfield(resultGUI,handles.SelectedDisplayOption)
    CubeNames = fieldnames(resultGUI);
    dose = resultGUI.(CubeNames{1,1});
else
    dose = resultGUI.(handles.SelectedDisplayOption);
end

%if function is called for the first time then set display parameters    
if isempty(handles.dispWindow{3,2})
    handles.dispWindow{3,1} = [min(dose(:)) max(dose(:))]; % set default dose range
    handles.dispWindow{3,2} = [min(dose(:)) max(dose(:))]; % set min max values  
end

minMaxRange = handles.dispWindow{3,1};
% if upper colorrange is defined then use it otherwise 120% iso dose 
 upperMargin = 1;
if abs((max(dose(:)) - handles.dispWindow{3,1}(1,2))) < 0.01  * max(dose(:))
    upperMargin = 1.2;
end

if (length(handles.IsoDose.Levels) == 1 && handles.IsoDose.Levels(1,1) == 0) || ~handles.IsoDose.NewIsoDoseFlag
    vLevels                  = [0.1:0.1:0.9 0.95:0.05:upperMargin];  
    referenceDose            = (minMaxRange(1,2))/(upperMargin);
    handles.IsoDose.Levels   = minMaxRange(1,1) + (referenceDose-minMaxRange(1,1)) * vLevels;
end
handles.IsoDose.Contours = matRad_computeIsoDoseContours(dose,handles.IsoDose.Levels);

    
   
%% CREATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% popup menu: machine
function popUpMachine_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% text box: max value
function txtMaxVal_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% text box: edit iso center
function editIsoCenter_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% text box: stratification levels
function editSequencingLevel_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slicerPrecision_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editBixelWidth_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editGantryAngle_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCouchAngle_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupRadMode_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editFraction_CreateFcn(hObject, ~, ~) 
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupPlane_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderSlice_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popupTypeOfPlot_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupDisplayOption_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderBeamSelection_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function listBoxCmd_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sliderOffset_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in legendTable.
function legendTable_Callback(hObject, eventdata, handles)
% hObject    handle to legendTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns legendTable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from legendTable
cst = evalin('base','cst');

idx    = get(hObject,'Value');
clr    = dec2hex(round(cst{idx,5}.visibleColor(:)*255),2)';
clr    = ['#';clr(:)]';

%Get the string entries
tmpString = get(handles.legendTable,'String');

if handles.VOIPlotFlag(idx)
    handles.VOIPlotFlag(idx) = false;
    cst{idx,5}.Visible = false;
    tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
elseif ~handles.VOIPlotFlag(idx)
    handles.VOIPlotFlag(idx) = true;
    cst{idx,5}.Visible = true;
    tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
end
set(handles.legendTable,'String',tmpString);

% update cst in workspace accordingly
assignin('base','cst',cst)

guidata(hObject, handles);
UpdatePlot(handles)

% --- Executes during object creation, after setting all properties.
function legendTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to legendTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in importDoseButton.
function importDoseButton_Callback(hObject,eventdata, handles)
% hObject    handle to importDoseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
extensions{1} = '*.nrrd';
[filenames,filepath,~] = uigetfile(extensions,'MultiSelect','on');

if ~iscell(filenames)
    tmp = filenames;
    filenames = cell(1);
    filenames{1} = tmp;
end

ct = evalin('base','ct');
resultGUI = evalin('base','resultGUI');

for filename = filenames
    [~,name,~] = fileparts(filename{1});
    [cube,~] = matRad_readCube(fullfile(filepath,filename{1}));
    if ~isequal(ct.cubeDim, size(cube))
        errordlg('Dimensions of the imported cube do not match with ct','Import failed!','modal');
        continue;
    end
    
    fieldname = ['import_' matlab.lang.makeValidName(name, 'ReplacementStyle','delete')];
    resultGUI.(fieldname) = cube;
end

assignin('base','resultGUI',resultGUI);
btnRefresh_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbutton_importFromBinary.
function pushbutton_importFromBinary_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_importFromBinary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % delete existing workspace - parse variables from base workspace
    set(handles.popupDisplayOption,'String','no option available');
    AllVarNames = evalin('base','who');
    RefVarNames = {'ct','cst','pln','stf','dij','resultGUI'};
    for i = 1:length(RefVarNames)  
        if sum(ismember(AllVarNames,RefVarNames{i}))>0
            evalin('base',['clear ', RefVarNames{i}]);
        end
    end
    handles.State = 0;

    %call the gui
    uiwait(matRad_importGUI);
    
    %Check if we have the variables in the workspace
    if evalin('base','exist(''cst'',''var'')') == 1 && evalin('base','exist(''ct'',''var'')') == 1
        cst = evalin('base','cst');
        ct = evalin('base','ct');
        %setCstTable(handles,cst);
        generateCstTable(hanles,cst);
        handles.TableChanged = false;
        set(handles.popupTypeOfPlot,'Value',1);
        
        % compute HU values
        if ~isfield(ct, 'cubeHU')
            ct = matRad_electronDensitiesToHU(ct);
            assignin('base','ct',ct);
        end
        if ~isfield(ct, 'cubeHU')
            handles.cubeHUavailable = false;
        else
            handles.cubeHUavailable = true;
        end
        
        % precompute contours 
        cst = precomputeContours(ct,cst);
    
        assignin('base','ct',ct);
        assignin('base','cst',cst);
        
        if evalin('base','exist(''pln'',''var'')')
            assignin('base','pln',pln);
            setPln(handles);
        else
            getPlnFromGUI(handles);
            setPln(handles);
        end
        handles.State = 1;
    end
    
    % set slice slider
    handles.plane = get(handles.popupPlane,'value');
    if handles.State >0
        set(handles.sliderSlice,'Min',1,'Max',ct.cubeDim(handles.plane),...
            'Value',round(ct.cubeDim(handles.plane)/2),...
            'SliderStep',[1/(ct.cubeDim(handles.plane)-1) 1/(ct.cubeDim(handles.plane)-1)]);
    end

    if handles.State > 0
        % define context menu for structures
        for i = 1:size(cst,1)
            if cst{i,5}.Visible
                handles.VOIPlotFlag(i) = true;
            else
                handles.VOIPlotFlag(i) = false;
            end
        end
    end
    
    handles.dispWindow = cell(3,2);
    handles.cBarChanged = true;
    
    UpdateState(handles);
    handles.rememberCurrAxes = false;
    UpdatePlot(handles);
    handles.rememberCurrAxes = true;
catch
   handles = showError(handles,'Binary Patient Import: Could not import data');
   UpdateState(handles);
end

guidata(hObject,handles);

% --- Executes on button press in radioBtnIsoCenter.
function radioBtnIsoCenter_Callback(hObject, eventdata, handles)
% hObject    handle to radioBtnIsoCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdatePlot(handles)
% Hint: get(hObject,'Value') returns toggle state of radioBtnIsoCenter

% --------------------------------------------------------------------
function uipushtool_screenshot_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_screenshot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
tmpFig = figure('position',[100 100 700 600],'Visible','off','name','Current View'); 
cBarHandle = findobj(handles.figure1,'Type','colorbar');
if ~isempty(cBarHandle)
    new_handle = copyobj([handles.axesFig cBarHandle],tmpFig);
else
    new_handle = copyobj(handles.axesFig,tmpFig);
end

oldPos = get(handles.axesFig,'Position');
set(new_handle(1),'units','normalized', 'Position',oldPos);

if ~isfield(handles,'lastStoragePath') || exist(handles.lastStoragePath,'dir') ~= 7
    handles.lastStoragePath = [];   
end

[filename, pathname] = uiputfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'; '*.fig','MATLAB figure file'},'Save current view',[handles.lastStoragePath 'screenshot.png']);

handles.lastStoragePath = pathname;

if ~isequal(filename,0) && ~isequal(pathname,0)
    set(gcf, 'pointer', 'watch');
    saveas(tmpFig,fullfile(pathname,filename));
    set(gcf, 'pointer', 'arrow');
    close(tmpFig);
    uiwait(msgbox('Current view has been succesfully saved!'));
else
    uiwait(msgbox('Aborted saving, showing figure instead!'));
    set(tmpFig,'Visible','on');
end

guidata(hObject,handles);


%% Callbacks & Functions for color setting
function UpdateColormapOptions(handles)

if isfield(handles,'colormapLocked') && handles.colormapLocked
    return;
end

selectionIndex = get(handles.popupmenu_chooseColorData,'Value');

cMapSelectionIndex = get(handles.popupmenu_chooseColormap,'Value');
cMapStrings = get(handles.popupmenu_chooseColormap,'String');

if selectionIndex > 1 
    set(handles.uitoggletool8,'State','on');
else
    set(handles.uitoggletool8,'State','off');
end

try 
    if selectionIndex == 2
        ct = evalin('base','ct');
        currentMap = handles.ctColorMap;
        window = handles.dispWindow{selectionIndex,1};
        if isfield(ct, 'cube')
            minMax = [min(ct.cube{1}(:)) max(ct.cube{1}(:))];
        else
            minMax = [min(ct.cubeHU{1}(:)) max(ct.cubeHU{1}(:))];
        end
        % adjust value for custom window to current
        handles.windowPresets(1).width = max(window) - min(window);
        handles.windowPresets(1).center = mean(window);
        % update full window information
        handles.windowPresets(2).width = minMax(2) - minMax(1);
        handles.windowPresets(2).center = mean(minMax);
    elseif selectionIndex == 3
        result = evalin('base','resultGUI');        
        dose = result.(handles.SelectedDisplayOption);
        currentMap = handles.doseColorMap;
        minMax = [min(dose(:)) max(dose(:))];
        window = handles.dispWindow{selectionIndex,1};
    else
        window = [0 1];
        minMax = window;
        currentMap = 'bone';
    end
catch
    window = [0 1];
    minMax = window;
    currentMap = 'bone';
end

valueRange = minMax(2) - minMax(1);

windowWidth = window(2) - window(1);
windowCenter = mean(window);

%This are some arbritrary settings to configure the sliders
sliderCenterMinMax = [minMax(1)-valueRange/2 minMax(2)+valueRange/2];
sliderWidthMinMax = [0 valueRange*2];

%if we have selected a value outside this range, we adapt the slider
%windows
if windowCenter < sliderCenterMinMax(1)
    sliderCenterMinMax(1) = windowCenter;
end
if windowCenter > sliderCenterMinMax(2)
    sliderCenterMinMax(2) = windowCenter;
end
if windowWidth < sliderWidthMinMax(1)
    sliderWidthMinMax(1) = windowWidth;
end
if windowCenter > sliderCenterMinMax(2)
    sliderWidthMinMax(2) = windowWidth;
end


set(handles.edit_windowCenter,'String',num2str(windowCenter,3));    
set(handles.edit_windowWidth,'String',num2str(windowWidth,3));
set(handles.edit_windowRange,'String',num2str(window,4));
set(handles.slider_windowCenter,'Min',sliderCenterMinMax(1),'Max',sliderCenterMinMax(2),'Value',windowCenter);
set(handles.slider_windowWidth,'Min',sliderWidthMinMax(1),'Max',sliderWidthMinMax(2),'Value',windowWidth);

cMapPopupIndex = find(strcmp(currentMap,cMapStrings));
set(handles.popupmenu_chooseColormap,'Value',cMapPopupIndex);

guidata(gcf,handles);

% --- Executes on selection change in popupmenu_chooseColorData.
function popupmenu_chooseColorData_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chooseColorData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_chooseColorData contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chooseColorData

%index = get(hObject,'Value') - 1;

handles.cBarChanged = true;

guidata(hObject,handles);
UpdatePlot(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_chooseColorData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chooseColorData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_windowCenter_Callback(hObject, eventdata, handles)
% hObject    handle to slider_windowCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newCenter      = get(hObject,'Value');
range          = get(handles.slider_windowWidth,'Value');
selectionIndex = get(handles.popupmenu_chooseColorData,'Value');

handles.dispWindow{selectionIndex,1}  = [newCenter-range/2 newCenter+range/2];
handles.cBarChanged = true;

guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function slider_windowCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_windowCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value',0.5);

% --- Executes on slider movement.
function slider_windowWidth_Callback(hObject, eventdata, handles)
% hObject    handle to slider_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newWidth = get(hObject,'Value');
center   = get(handles.slider_windowCenter,'Value');
selectionIndex                        = get(handles.popupmenu_chooseColorData,'Value');
handles.dispWindow{selectionIndex,1}  = [center-newWidth/2 center+newWidth/2];
handles.cBarChanged = true;

guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function slider_windowWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

set(hObject,'Value',1.0);


% --- Executes on selection change in popupmenu_chooseColormap.
function popupmenu_chooseColormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_chooseColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_chooseColormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_chooseColormap

index = get(hObject,'Value');
strings = get(hObject,'String');

selectionIndex = get(handles.popupmenu_chooseColorData,'Value');

switch selectionIndex 
    case 2
        handles.ctColorMap = strings{index};
    case 3
        handles.doseColorMap = strings{index};
    otherwise
end

handles.cBarChanged = true;

guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_chooseColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_chooseColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_windowRange_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowRange as text
%        str2double(get(hObject,'String')) returns contents of edit_windowRange as a double

selectionIndex = get(handles.popupmenu_chooseColorData,'Value');
vRange         = str2num(get(hObject,'String'));
% matlab adds a zero in the beginning when text field is changed
if numel(vRange) == 3
    vRange = vRange(vRange~=0);
end

handles.dispWindow{selectionIndex,1} = sort(vRange);

handles.cBarChanged = true;

    % compute new iso dose lines
if selectionIndex > 2 
    guidata(hObject,handles);
    handles = updateIsoDoseLineCache(handles);
end

guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function edit_windowRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0 1');


function edit_windowCenter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowCenter as text
%        str2double(get(hObject,'String')) returns contents of edit_windowCenter as a double

newCenter           = str2double(get(hObject,'String'));
width               = get(handles.slider_windowWidth,'Value');
selectionIndex      = get(handles.popupmenu_chooseColorData,'Value');
handles.dispWindow{selectionIndex,1}  = [newCenter-width/2 newCenter+width/2];
handles.cBarChanged = true;
guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function edit_windowCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_windowWidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowWidth as text
%        str2double(get(hObject,'String')) returns contents of edit_windowWidth as a double

newWidth            = str2double(get(hObject,'String'));
center              = get(handles.slider_windowCenter,'Value');
selectionIndex      = get(handles.popupmenu_chooseColorData,'Value');
handles.dispWindow{selectionIndex,1}  = [center-newWidth/2 center+newWidth/2];
handles.cBarChanged = true;
guidata(hObject,handles);
UpdatePlot(handles);


% --- Executes during object creation, after setting all properties.
function edit_windowWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uitoggletool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if on or off
val = strcmp(get(hObject,'State'),'on');

%Now we have to apply the new selection to our colormap options panel
if ~val
    newSelection = 1;
else
    %Chooses the selection from the highest state
    selections = get(handles.popupmenu_chooseColorData,'String');
    newSelection = numel(selections);
end    
set(handles.popupmenu_chooseColorData,'Value',newSelection);

handles.cBarChanged = true;
guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes on slider movement.
function sliderOpacity_Callback(hObject, eventdata, handles)
% hObject    handle to sliderOpacity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.doseOpacity = get(hObject,'Value');
guidata(hObject,handles);
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function sliderOpacity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderOpacity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%% Data Cursors
function cursorText = dataCursorUpdateFunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

target = findall(0,'Name','matRadGUI');

% Get GUI data (maybe there is another way?)
handles = guidata(target);

% position of the data point to label
pos = get(event_obj,'Position');

%Different behavior for image and profile plot
if get(handles.popupTypeOfPlot,'Value')==1 %Image view
    cursorText = cell(0,1);
    try   
        if handles.State >= 1
            plane = get(handles.popupPlane,'Value');
            slice = round(get(handles.sliderSlice,'Value'));
            
            %Get the CT values
            ct  = evalin('base','ct');
            
            %We differentiate between pos and ix, since the user may put
            %the datatip on an isoline which returns a continous position
            cubePos = zeros(1,3);
            cubePos(plane) = slice;
            cubePos(1:end ~= plane) = fliplr(pos);            
            cubeIx = round(cubePos);
            
            %Here comes the index permutation stuff
            %Cube Index
            cursorText{end+1,1} = ['Cube Index: ' mat2str(cubeIx)];
            %Space Coordinates
            coords = zeros(1,3);
            coords(1) = cubePos(2)*ct.resolution.y;
            coords(2) = cubePos(1)*ct.resolution.x;
            coords(3) = cubePos(3)*ct.resolution.z;            
            cursorText{end+1,1} = ['Space Coordinates: ' mat2str(coords,5) ' mm'];
            
            ctVal = ct.cubeHU{1}(cubeIx(1),cubeIx(2),cubeIx(3));
            cursorText{end+1,1} = ['HU Value: ' num2str(ctVal,3)];
        end
        
        %Add dose information if available
        if handles.State == 3
            %get result structure
            result = evalin('base','resultGUI');
            
            %Get all result names from popup
            resultNames = get(handles.popupDisplayOption,'String');
            
            %Display all values of fields found in the resultGUI struct
            for runResult = 1:numel(resultNames)               
                name = resultNames{runResult};
                if isfield(result,name)
                    field = result.(name);
                    val = field(cubeIx(1),cubeIx(2),cubeIx(3));
                    cursorText{end+1,1} = [name ': ' num2str(val,3)];
                end
            end      
        end
    catch
        cursorText{end+1,1} = 'Error while retreiving Data!';
    end    
else %Profile view
    cursorText = cell(2,1);
    cursorText{1} = ['Radiological Depth: ' num2str(pos(1),3) ' mm'];
    cursorText{2} = [get(target,'DisplayName') ': ' num2str(pos(2),3)];
end



% --- Executes on selection change in popMenuBioOpt.
function popMenuBioOpt_Callback(hObject, ~, handles)
pln = evalin('base','pln');
contentBioOpt = get(handles.popMenuBioOpt,'String');
NewBioOptimization = contentBioOpt(get(handles.popMenuBioOpt,'Value'),:);

if handles.State > 0
    if (strcmp(pln.propOpt.bioOptimization,'LEMIV_effect') && strcmp(NewBioOptimization,'LEMIV_RBExD')) ||...
       (strcmp(pln.propOpt.bioOptimization,'LEMIV_RBExD') && strcmp(NewBioOptimization,'LEMIV_effect')) 
       % do nothing - re-optimization is still possible
    elseif ((strcmp(pln.propOpt.bioOptimization,'const_RBE') && strcmp(NewBioOptimization,'none')) ||...
           (strcmp(pln.propOpt.bioOptimization,'none') && strcmp(NewBioOptimization,'const_RBE'))) && isequal(pln.radiationMode,'protons')
       % do nothing - re-optimization is still possible
    else
        handles.State = 1;
    end
end
getPlnFromGUI(handles);

UpdateState(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popMenuBioOpt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn3Dview.
function btn3Dview_Callback(hObject, eventdata, handles)
% hObject    handle to btn3Dview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'axesFig3D') || ~isfield(handles,'axesFig3D') || ~isgraphics(handles.axesFig3D)
    handles.fig3D = figure('Name','matRad 3D View');
    handles.axesFig3D = axes('Parent',handles.fig3D);
    view(handles.axesFig3D,3);
end
%end

UpdatePlot(handles);

guidata(hObject,handles);

% --- Executes on button press in radiobtnCT.
function radiobtnCT_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnCT
UpdatePlot(handles)

% --- Executes on button press in radiobtnPlan.
function radiobtnPlan_Callback(hObject, eventdata, handles)
% hObject    handle to radiobtnPlan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobtnPlan
UpdatePlot(handles)


% --- Executes on selection change in popupmenu_windowPreset.
function popupmenu_windowPreset_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_windowPreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_windowPreset contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_windowPreset

selectionIndexCube      = 2; % working on ct only
selectionIndexWindow    = get(handles.popupmenu_windowPreset,'Value');
newCenter               = handles.windowPresets(selectionIndexWindow).center;
newWidth                = handles.windowPresets(selectionIndexWindow).width;

handles.dispWindow{selectionIndexCube,1}  = [newCenter - newWidth/2 newCenter + newWidth/2];
handles.cBarChanged = true;
guidata(hObject,handles);
UpdatePlot(handles);
UpdateColormapOptions(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_windowPreset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_windowPreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% setup ct window list
% data and values from CERR https://github.com/adityaapte/CERR
windowNames = {'Custom','Full','Abd/Med', 'Head', 'Liver', 'Lung', 'Spine', 'Vrt/Bone'};
windowCenter = {NaN, NaN, -10, 45, 80, -500, 30, 400};
windowWidth = {NaN, NaN, 330, 125, 305, 1500, 300, 1500};
windowPresets = cell2struct([windowNames', windowCenter', windowWidth'], {'name', 'center', 'width'},2);


handles.windowPresets = windowPresets;

selectionList = {windowPresets(:).name};
set(hObject,'String',selectionList(:));
set(hObject,'Value',1);


guidata(hObject,handles);

% --- Executes on button press in checkbox_lockColormap.
function checkbox_lockColormap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_lockColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_lockColormap
handles.colormapLocked = get(hObject,'Value');

if handles.colormapLocked
    state = 'Off'; %'Inactive';
else
    state = 'On';
end

set(handles.popupmenu_chooseColorData,'Enable',state);
set(handles.popupmenu_windowPreset,'Enable',state);
set(handles.slider_windowWidth,'Enable',state);
set(handles.slider_windowCenter,'Enable',state);
set(handles.edit_windowWidth,'Enable',state);
set(handles.edit_windowCenter,'Enable',state);
set(handles.edit_windowRange,'Enable',state);
set(handles.popupmenu_chooseColormap,'Enable',state);


guidata(hObject,handles);


function cst = updateStructureTable(handles,cst)
colorAssigned = true;

% check whether all structures have an assigned color
for i = 1:size(cst,1)
    if ~isfield(cst{i,5},'visibleColor')
        colorAssigned = false;
        break;
    elseif isempty(cst{i,5}.visibleColor)
        colorAssigned = false;
        break;
    end
end

% assign color if color assignment is not already present or inconsistent
if colorAssigned == false
  m         = 64;
  colorStep = ceil(m/size(cst,1));
  colors    = colorcube(colorStep*size(cst,1));
  % spread individual VOI colors in the colorcube color palette
  colors    = colors(1:colorStep:end,:);
  
  for i = 1:size(cst,1)
    cst{i,5}.visibleColor = colors(i,:);
  end
end

for s = 1:size(cst,1)
    handles.VOIPlotFlag(s) = cst{s,5}.Visible;
    clr = dec2hex(round(cst{s,5}.visibleColor(:)*255),2)';
    clr = ['#';clr(:)]';
    if handles.VOIPlotFlag(s)
        tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
    else
        tmpString{s} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{s,2},'</TD></TR> </table></html>'];
    end
end
set(handles.legendTable,'String',tmpString);

%-- generates the CST table
function cst = generateCstTable(handles,cst)

cst = updateStructureTable(handles,cst);

cstPanel = handles.uipanel3;

cstPanelPos = getpixelposition(cstPanel);

%Parameters for line height
objHeight = 22;
lineHeight = 25; %Height of a table line
fieldSep = 2; %Separation between fields horizontally
yTopSep = 40; %Separation of the first line from the top
tableViewHeight = cstPanelPos(4) - yTopSep; %Full height of the view

%Widths of the fields
buttonW = objHeight;
nameW = 90;
typeW = 70;
opW = objHeight;
functionW = 120;
penaltyW = 40;
paramTitleW = 120;
paramW = 30;


%Scrollbar
cstVertTableScroll = findobj(cstPanel.Children,'Style','slider');
if isempty(cstVertTableScroll)
    sliderPos = 0;
else
    sliderPos = cstVertTableScroll.Max - cstVertTableScroll.Value;
end
%disp(num2str(sliderPos));
ypos = @(c) tableViewHeight - c*lineHeight + sliderPos;

delete(cstPanel.Children);

%Creates a dummy axis to allow for the use of textboxes instead of uicontrol to be able to use the (la)tex interpreter
tmpAxes = axes('Parent',cstPanel,'units','normalized','position',[0 0 1 1],'visible','off');

organTypes = {'OAR', 'TARGET'};

%columnname = {'VOI name','VOI type','priority','obj. / const.'};%,'penalty','dose', 'EUD','volume','robustness'};
           
%Get all Classes & classNames                   
mpkgObjectives = meta.package.fromName('DoseObjectives');
mpkgConstraints = meta.package.fromName('DoseConstraints');
classList = [mpkgObjectives.ClassList; mpkgConstraints.ClassList];

classList = classList(not([classList.Abstract]));

%Now get the "internal" name from the properties
classNames = cell(2,numel(classList));
for clIx = 1:numel(classList)
    cl = classList(clIx);
    pList = cl.PropertyList; %Get List of all properties
    pNameIx = arrayfun(@(p) strcmp(p.Name,'name'),pList); %get index of the "name" property
    p = pList(pNameIx); %select name property
    pName = p.DefaultValue; % get value / name
    classNames(:,clIx) = {cl.Name; pName}; %Store class name and display name
end

% Collect Class-File & Display Names
%classNames = {classList.Name; p.DefaultValue};

%columnformat = {cst(:,2)',{'OAR','TARGET'},'numeric',...
%       AllObjectiveFunction,...
%       'numeric','numeric','numeric','numeric',{'none','WC','prob'}};
   
numOfObjectives = 0;
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        numOfObjectives = numOfObjectives + numel(cst{i,6});
    end
end

cnt = 0;

newline = '\n';

%Setup Headlines
xPos = 5;
h = uicontrol(cstPanel,'Style','text','String','+/-','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Remove or add Constraint or Objective');
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','VOI name','Position',[xPos ypos(cnt) nameW objHeight],'TooltipString','Name of the structure with objective/constraint');
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','VOI type','Position',[xPos ypos(cnt) typeW objHeight],'TooltipString','Segmentation Classification');
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','OP','Position',[xPos ypos(cnt) opW objHeight],'TooltipString',['Overlap Priority' char(10) '(Smaller number overlaps higher number)']);
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','Function','Position',[xPos ypos(cnt) functionW objHeight],'TooltipString','Objective/Constraint function type');
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','p','Position',[xPos ypos(cnt) penaltyW objHeight],'TooltipString','Optimization penalty');
xPos = xPos + h.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','text','String','| Parameters','Position',[xPos ypos(cnt) paramTitleW objHeight],'TooltipString','List of parameters','HorizontalAlignment','left');
xPos = xPos + h.Position(3) + fieldSep;
cnt = cnt + 1;

%Create Objectives / Constraints controls
for i = 1:size(cst,1)   
   if strcmp(cst(i,3),'IGNORED')~=1
      for j=1:numel(cst{i,6})
           
           obj = cst{i,6}{j};
           
           %Convert to class if not
           if ~isa(obj,'matRad_DoseOptimizationFunction')
                try
                    obj = matRad_DoseOptimizationFunction.createInstanceFromStruct(obj);
                catch
                    warning('Objective/Constraint not valid!')
                    continue;
                end
           end
               
           %VOI
           %data{Counter,1}  = cst{i,2};
           %ypos = cstPanelPos(4) - (yTopSep + cnt*lineHeight);
           xPos = 5;
           
           %h = uicontrol(cstPanel,'Style','popupmenu','String',cst(:,2)','Position',[xPos ypos 100 objHeight]);
           %h.Value = i;
           h = uicontrol(cstPanel,'Style','pushbutton','String','-','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Remove Objective/Constraint','Callback',{@btObjRemove_Callback,handles},...
               'UserData',[i,j]);
           xPos = xPos + h.Position(3) + fieldSep;
           h = uicontrol(cstPanel','Style','edit','String',cst{i,2},'Position',[xPos ypos(cnt) nameW objHeight],'TooltipString','Name',...
               'Enable','inactive',... %Disable editing of name atm
               'UserData',[i,2],'Callback',{@editCstParams_Callback,handles}); %Callback added, however, editing is disabled atm
           xPos = xPos + h.Position(3) + fieldSep;
           h = uicontrol(cstPanel,'Style','popupmenu','String',organTypes','Value',find(strcmp(cst{i,3},organTypes)),'Position',[xPos ypos(cnt) typeW objHeight],'TooltipString','Segmentation Classification',...
               'UserData',[i,3],'Callback',{@editCstParams_Callback,handles});
           xPos = xPos + h.Position(3) + fieldSep;
           h = uicontrol(cstPanel,'Style','edit','String',num2str(cst{i,5}.Priority),'Position',[xPos ypos(cnt) opW objHeight],'TooltipString',['Overlap Priority' newline '(Smaller number overlaps higher number)'],...
               'UserData',[i,5],'Callback',{@editCstParams_Callback,handles});
           xPos = xPos + h.Position(3) + fieldSep;
           
           h = uicontrol(cstPanel,'Style','popupmenu','String',classNames(2,:)','Value',find(strcmp(obj.name,classNames(2,:))),'Position',[xPos ypos(cnt) functionW objHeight],'TooltipString','Select Objective/Constraint',...
               'UserData',{[i,j],classNames(1,:)},'Callback',{@changeObjFunction_Callback,handles});
           xPos = xPos + h.Position(3) + fieldSep;
           
           %Check if we have an objective to display penalty
           if isa(obj,'DoseObjectives.matRad_DoseObjective')
               h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.penalty),'Position',[xPos ypos(cnt) penaltyW objHeight],'TooltipString','Objective Penalty','UserData',[i,j,0],'Callback',{@editObjParam_Callback,handles});
           else
               h = uicontrol(cstPanel,'Style','edit','String','----','Position',[xPos ypos(cnt) penaltyW objHeight],'Enable','off');
           end
           xPos = xPos + h.Position(3) + fieldSep;
           
           for p = 1:numel(obj.parameterNames)
              %h = uicontrol(cstPanel,'Style','edit','String',obj.parameters{1,p},'Position',[xPos ypos(cnt) 100 objHeight],'Enable','inactive');
              %xPos = xPos + h.Position(3) + fieldSep;
              h = text('Parent',tmpAxes,'String',['| ' obj.parameterNames{p} ':'],'VerticalAlignment','middle','units','pix','Position',[xPos ypos(cnt)+lineHeight/2],'Interpreter','tex','FontWeight','normal',...
                  'FontSize',cstPanel.FontSize,'FontName',cstPanel.FontName,'FontUnits',cstPanel.FontUnits,'FontWeight','normal');%[xPos ypos(cnt) 100 objHeight]);
              xPos = xPos + h.Extent(3) + fieldSep;
              %h = annotation(cstPanel,'textbox','String',obj.parameters{1,p},'Units','pix','Position', [xPos ypos(cnt) 100 objHeight],'Interpreter','Tex');
              
              %Check if we have a cell and therefore a parameter list
              if iscell(obj.parameterTypes{p})                  
                  h = uicontrol(cstPanel,'Style','popupmenu','String',obj.parameterTypes{p}','Value',obj.parameters{p},'TooltipString',obj.parameterNames{p},'Position',[xPos ypos(cnt) paramW*2 objHeight],'UserData',[i,j,p],'Callback',{@editObjParam_Callback,handles});
              else
                h = uicontrol(cstPanel,'Style','edit','String',num2str(obj.parameters{p}),'TooltipString',obj.parameterNames{p},'Position',[xPos ypos(cnt) paramW objHeight],'UserData',[i,j,p],'Callback',{@editObjParam_Callback,handles});
              end
              
              xPos = xPos + h.Position(3) + fieldSep;
           end

           cnt = cnt +1;
       end
   end   
end
xPos = 5;
hAdd = uicontrol(cstPanel,'Style','pushbutton','String','+','Position',[xPos ypos(cnt) buttonW objHeight],'TooltipString','Add Objective/Constraint','Callback',{@btObjAdd_Callback,handles});
xPos = xPos + hAdd.Position(3) + fieldSep;
h = uicontrol(cstPanel,'Style','popupmenu','String',cst(:,2)','Position',[xPos ypos(cnt) nameW objHeight]);
hAdd.UserData = h;

%Calculate Scrollbar
lastPos = ypos(cnt);
firstPos = ypos(0);
tableHeight = abs(firstPos - lastPos);

exceedFac = tableHeight / tableViewHeight;
if exceedFac > 1
    sliderFac = exceedFac - 1; 
    uicontrol(cstPanel,'Style','slider','Units','normalized','Position',[0.975 0 0.025 1],'Min',0,'Max',ceil(sliderFac)*tableViewHeight,'SliderStep',[lineHeight tableViewHeight] ./ (ceil(sliderFac)*tableViewHeight),'Value',ceil(sliderFac)*tableViewHeight - sliderPos,'Callback',{@cstTableSlider_Callback,handles});
end




%set(handles.uiTable,'ColumnName',columnname);
%set(handles.uiTable,'ColumnFormat',columnformat);
%set(handles.uiTable,'ColumnEditable',[true true true true true true true true true true]);
%set(handles.uiTable,'Data',data);
   

% --- Executes when uipanel3 is resized.
function uipanel3_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    generateCstTable(handles,evalin('base','cst'));
catch
end

function btObjAdd_Callback(hObject, ~, handles)
% hObject    handle to btnuiTableAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupHandle = hObject.UserData;
cstIndex = popupHandle.Value;

cst = evalin('base','cst');
%Add Standard Objective
if strcmp(cst{cstIndex,3},'TARGET')
    cst{cstIndex,6}{end+1} = struct(DoseObjectives.matRad_SquaredDeviation);
else
    cst{cstIndex,6}{end+1} = struct(DoseObjectives.matRad_SquaredOverdosing);
end

assignin('base','cst',cst);

%set(handles.uiTable,'data',data);

%handles.State=1;
%guidata(hObject,handles);
%UpdateState(handles);

generateCstTable(handles,cst);

function btObjRemove_Callback(hObject, ~, handles)
% hObject    handle to btnuiTableAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ix = hObject.UserData;

cst = evalin('base','cst');
%Add Standard Objective

cst{ix(1),6}(ix(2)) = [];

assignin('base','cst',cst);

%set(handles.uiTable,'data',data);

%handles.State=1;
%guidata(hObject,handles);
%UpdateState(handles);

generateCstTable(handles,cst);

function editObjParam_Callback(hObject, ~, handles)
% hObject    handle to btnuiTableAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ix = hObject.UserData;

cst = evalin('base','cst');
%Add Standard Objective

%if the third index is 0 we changed the penalty
%if we have a popupmenu selection we use value
%otherwise we use the edit string

%{
%obj = cst{ix(1),6}{ix(2)};
%Convert to class if not
if ~isa(obj,'matRad_DoseOptimizationFunction')
    try
        eval([obj.className '(obj)']);
    catch
        warning('Objective/Constraint not valid!')
        return;
    end
end
%}

if ix(3) == 0
    cst{ix(1),6}{ix(2)}.penalty = str2double(hObject.String);
elseif isequal(hObject.Style,'popupmenu')
    cst{ix(1),6}{ix(2)}.parameters{ix(3)} = hObject.Value;
else
    cst{ix(1),6}{ix(2)}.parameters{ix(3)} = str2double(hObject.String);
end
    
assignin('base','cst',cst);

%set(handles.uiTable,'data',data);

%handles.State=1;
%guidata(hObject,handles);
%UpdateState(handles);

generateCstTable(handles,cst);

function changeObjFunction_Callback(hObject, ~, handles)
% hObject    handle to btnuiTableAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = hObject.UserData;
ix = data{1};
classNames = data{2};
classToCreate = classNames{hObject.Value};

cst = evalin('base','cst');
%Add Standard Objective

%We just check if the user really wanted to change the objective to be
%user-friendly
currentObj = cst{ix(1),6}{ix(2)};
currentClass = class(currentObj);
if ~strcmp(currentClass,classToCreate)    
    newObj = eval(classToCreate);
    
    % Only if we have a penalty value for optimization, apply the new one
    % Maybe this check should be more exact?        
    if (isfield(currentObj,'penalty') || isprop(currentObj,'penalty')) && isprop(newObj,'penalty')
        newObj.penalty = currentObj.penalty;
    end
    
    cst{ix(1),6}{ix(2)} = struct(newObj);
    
    assignin('base','cst',cst);
    
    %set(handles.uiTable,'data',data);
    
    %handles.State=1;
    %guidata(hObject,handles);
    %UpdateState(handles);
    
    generateCstTable(handles,cst);
end

function editCstParams_Callback(hObject,~,handles)
data = hObject.UserData;
ix = data(1);
col = data(2);

cst = evalin('base','cst');

switch col
    case 2
        cst{ix,col} = hObject.String;
    case 3
        cst{ix,col} = hObject.String{hObject.Value};
    case 5
        cst{ix,col}.Priority = uint32(str2double(hObject.String));
    otherwise
        warning('Wrong column assignment in GUI based cst setting');
end

assignin('base','cst',cst);

generateCstTable(handles,cst);

    
% --- Executes on button press in radiobutton3Dconf.
function radiobutton3Dconf_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3Dconf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3Dconf

function x = matRad_checkForConnectedBixelRows(stf)

x = true;

for i = 1:size(stf,2)
    
    bixelPos = reshape([stf(i).ray.rayPos_bev],3,[]);
   
    rowCoords = unique(bixelPos(3,:));
    
    for j = 1:numel(rowCoords)
        
        increments = diff(bixelPos(1,rowCoords(j) == bixelPos(3,:)));
        
        % if we find one not connected row -> return false
        if numel(unique(increments)) > 1
            x = false;
            return;
        end
    end
    
end

% --- Executes on slider movement.
function cstTableSlider_Callback(hObject, eventdata, handles)
% hObject    handle to cstTableSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
generateCstTable(handles,evalin('base','cst'));

% --- Executes during object creation, after setting all properties.
function cstTableSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cstTableSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Creates and returns a handle to the GUI figure. 
function h1 = matRadGUI_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end
load matRadGUI_export.mat


appdata = [];
appdata.GUIDEOptions = mat{1};
appdata.lastValidTag = 'figure1';
appdata.UsedByGUIData_m = struct(...
    'windowPresets', struct(...
    'name', 'Custom', ...
    'center', NaN, ...
    'width', NaN));
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Units','characters',...
'Position',[138.4 -7.38461538461539 273.4 59.5384615384615],...
'Visible','on',...
'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
'CloseRequestFcn',@(hObject,eventdata)matRadGUI_export('figure1_CloseRequestFcn',hObject,eventdata,guidata(hObject)),...
'IntegerHandle','off',...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'MenuBar','none',...
'Name','matRadGUI',...
'NumberTitle','off',...
'Tag','figure1',...
'UserData',[],...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.99999864 29.69999902],...
'PaperType',get(0,'defaultfigurePaperType'),...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'ScreenPixelsPerInchMode','manual',...
'HandleVisibility','callback',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'axesLogo';

h2 = axes(...
'Parent',h1,...
'FontUnits',get(0,'defaultaxesFontUnits'),...
'Units',get(0,'defaultaxesUnits'),...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'CameraTarget',[0.5 0.5 0.5],...
'CameraTargetMode',get(0,'defaultaxesCameraTargetMode'),...
'CameraViewAngle',6.60861036031192,...
'CameraViewAngleMode',get(0,'defaultaxesCameraViewAngleMode'),...
'PlotBoxAspectRatio',[1 0.288951841359773 0.288951841359773],...
'PlotBoxAspectRatioMode',get(0,'defaultaxesPlotBoxAspectRatioMode'),...
'FontName','CMU Serif',...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'ColormapMode',get(0,'defaultaxesColormapMode'),...
'Alphamap',[0 0.0159 0.0317 0.0476 0.0635 0.0794 0.0952 0.1111 0.127 0.1429 0.1587 0.1746 0.1905 0.2063 0.2222 0.2381 0.254 0.2698 0.2857 0.3016 0.3175 0.3333 0.3492 0.3651 0.381 0.3968 0.4127 0.4286 0.4444 0.4603 0.4762 0.4921 0.5079 0.5238 0.5397 0.5556 0.5714 0.5873 0.6032 0.619 0.6349 0.6508 0.6667 0.6825 0.6984 0.7143 0.7302 0.746 0.7619 0.7778 0.7937 0.8095 0.8254 0.8413 0.8571 0.873 0.8889 0.9048 0.9206 0.9365 0.9524 0.9683 0.9841 1],...
'AlphamapMode',get(0,'defaultaxesAlphamapMode'),...
'XTick',[0 0.2 0.4 0.6 0.8 1],...
'XTickMode',get(0,'defaultaxesXTickMode'),...
'XTickLabel',{  '0'; '0.2'; '0.4'; '0.6'; '0.8'; '1' },...
'XTickLabelMode',get(0,'defaultaxesXTickLabelMode'),...
'YTick',[0 0.5 1],...
'YTickMode',get(0,'defaultaxesYTickMode'),...
'YTickLabel',{  '0'; '0.5'; '1' },...
'YTickLabelMode',get(0,'defaultaxesYTickLabelMode'),...
'Color',get(0,'defaultaxesColor'),...
'CameraMode',get(0,'defaultaxesCameraMode'),...
'DataSpaceMode',get(0,'defaultaxesDataSpaceMode'),...
'ColorSpaceMode',get(0,'defaultaxesColorSpaceMode'),...
'DecorationContainerMode',get(0,'defaultaxesDecorationContainerMode'),...
'ChildContainerMode',get(0,'defaultaxesChildContainerMode'),...
'XRulerMode',get(0,'defaultaxesXRulerMode'),...
'YRulerMode',get(0,'defaultaxesYRulerMode'),...
'ZRulerMode',get(0,'defaultaxesZRulerMode'),...
'AmbientLightSourceMode',get(0,'defaultaxesAmbientLightSourceMode'),...
'Tag','axesLogo',...
'Position',[0.444874274661509 0.893725992317542 0.184397163120567 0.0998719590268886],...
'ActivePositionProperty','position',...
'ActivePositionPropertyMode',get(0,'defaultaxesActivePositionPropertyMode'),...
'LooseInset',[0.182759687929063 0.112926163636008 0.133555156563546 0.0769951115700055],...
'FontSize',9.63,...
'SortMethod','childorder',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h3 = get(h2,'title');

set(h3,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0 0 0],...
'ColorMode','auto',...
'Position',[0.500002778623327 1.02596323529412 0.5],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','bold',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','Axes Title',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h4 = get(h2,'xlabel');

set(h4,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0.500000476837158 -0.277303926477245 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','top',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h5 = get(h2,'ylabel');

set(h5,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[-0.0943626077857305 0.500000476837158 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',90,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h6 = get(h2,'zlabel');

set(h6,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0 0 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','left',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','middle',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','off',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

appdata = [];
appdata.lastValidTag = 'axesDKFZ';

h7 = axes(...
'Parent',h1,...
'FontUnits',get(0,'defaultaxesFontUnits'),...
'Units',get(0,'defaultaxesUnits'),...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'CameraTarget',[0.5 0.5 0.5],...
'CameraTargetMode',get(0,'defaultaxesCameraTargetMode'),...
'CameraViewAngle',6.60861036031192,...
'CameraViewAngleMode',get(0,'defaultaxesCameraViewAngleMode'),...
'PlotBoxAspectRatio',[1 0.145161290322581 0.145161290322581],...
'PlotBoxAspectRatioMode',get(0,'defaultaxesPlotBoxAspectRatioMode'),...
'FontName','CMU Serif',...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'ColormapMode',get(0,'defaultaxesColormapMode'),...
'Alphamap',[0 0.0159 0.0317 0.0476 0.0635 0.0794 0.0952 0.1111 0.127 0.1429 0.1587 0.1746 0.1905 0.2063 0.2222 0.2381 0.254 0.2698 0.2857 0.3016 0.3175 0.3333 0.3492 0.3651 0.381 0.3968 0.4127 0.4286 0.4444 0.4603 0.4762 0.4921 0.5079 0.5238 0.5397 0.5556 0.5714 0.5873 0.6032 0.619 0.6349 0.6508 0.6667 0.6825 0.6984 0.7143 0.7302 0.746 0.7619 0.7778 0.7937 0.8095 0.8254 0.8413 0.8571 0.873 0.8889 0.9048 0.9206 0.9365 0.9524 0.9683 0.9841 1],...
'AlphamapMode',get(0,'defaultaxesAlphamapMode'),...
'XTick',[0 0.2 0.4 0.6 0.8 1],...
'XTickMode',get(0,'defaultaxesXTickMode'),...
'XTickLabel',{  '0'; '0.2'; '0.4'; '0.6'; '0.8'; '1' },...
'XTickLabelMode',get(0,'defaultaxesXTickLabelMode'),...
'YTick',[0 0.5 1],...
'YTickMode',get(0,'defaultaxesYTickMode'),...
'YTickLabel',{  '0'; '0.5'; '1' },...
'YTickLabelMode',get(0,'defaultaxesYTickLabelMode'),...
'Color',get(0,'defaultaxesColor'),...
'CameraMode',get(0,'defaultaxesCameraMode'),...
'DataSpaceMode',get(0,'defaultaxesDataSpaceMode'),...
'ColorSpaceMode',get(0,'defaultaxesColorSpaceMode'),...
'DecorationContainerMode',get(0,'defaultaxesDecorationContainerMode'),...
'ChildContainerMode',get(0,'defaultaxesChildContainerMode'),...
'XRulerMode',get(0,'defaultaxesXRulerMode'),...
'YRulerMode',get(0,'defaultaxesYRulerMode'),...
'ZRulerMode',get(0,'defaultaxesZRulerMode'),...
'AmbientLightSourceMode',get(0,'defaultaxesAmbientLightSourceMode'),...
'Tag','axesDKFZ',...
'Position',[0.652482269503546 0.903969270166453 0.259187620889749 0.0717029449423816],...
'ActivePositionProperty','position',...
'ActivePositionPropertyMode',get(0,'defaultaxesActivePositionPropertyMode'),...
'LooseInset',[0.13 0.11 0.0950000000000001 0.0749999999999999],...
'FontSize',9.63,...
'SortMethod','childorder',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h8 = get(h7,'title');

set(h8,...
'Parent',h7,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0 0 0],...
'ColorMode','auto',...
'Position',[0.500002880250254 1.03678125 0.500000000000007],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','bold',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','Axes Title',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h9 = get(h7,'xlabel');

set(h9,...
'Parent',h7,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0.500000476837158 -0.392847229176095 7.105427357601e-15],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','top',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h10 = get(h7,'ylabel');

set(h10,...
'Parent',h7,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[-0.0671572591700866 0.50000047683716 7.105427357601e-15],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',90,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h11 = get(h7,'zlabel');

set(h11,...
'Parent',h7,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0 0 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','left',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','middle',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','off',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

appdata = [];
appdata.lastValidTag = 'uipanel1';

h12 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'ShadowColor',get(0,'defaultuipanelShadowColor'),...
'Title','Plan',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel1',...
'UserData',[],...
'Clipping','off',...
'Position',[0.00451321727917473 0.527528809218956 0.430689877498388 0.272727272727273],...
'FontName','Helvetica',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtBixelWidth';

h13 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','bixel width in [mm]',...
'Style','text',...
'Position',[0.0347043701799486 0.859315589353612 0.17866323907455 0.0950570342205324],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtBixelWidth',...
'UserData',[],...
'FontName','Helvetica',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editBixelWidth';

h14 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','5',...
'Style','edit',...
'Position',[0.219794344473008 0.889733840304182 0.161953727506427 0.0836501901140684],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editBixelWidth_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editBixelWidth_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editBixelWidth',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtGantryAngle';

h15 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Gantry Angle in °',...
'Style','text',...
'Position',[0.032133676092545 0.752851711026616 0.176092544987147 0.0950570342205324],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtGantryAngle',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editGantryAngle';

h16 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','0',...
'Style','edit',...
'Position',[0.219794344473008 0.779467680608365 0.161953727506427 0.0836501901140684],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editGantryAngle_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editGantryAngle_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editGantryAngle',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtCouchAngle';

h17 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Couch Angle in °',...
'Style','text',...
'Position',[0.0347043701799486 0.64638783269962 0.173521850899743 0.0950570342205324],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtCouchAngle',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editCouchAngle';

h18 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','0',...
'Style','edit',...
'Position',[0.219794344473008 0.669201520912547 0.161953727506427 0.0836501901140685],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editCouchAngle_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editCouchAngle_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editCouchAngle',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'popupRadMode';

h19 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'photons'; 'protons'; 'carbon' },...
'Style','popupmenu',...
'Value',1,...
'Position',[0.219794344473008 0.52851711026616 0.161953727506427 0.114068441064639],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupRadMode_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupRadMode_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupRadMode',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtRadMode';

h20 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Radiation Mode',...
'Style','text',...
'Position',[0.051413881748072 0.539923954372624 0.146529562982005 0.0950570342205324],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtRadMode',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtNumOfFractions';

h21 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','# Fractions',...
'Style','text',...
'Position',[0.051413881748072 0.209125475285171 0.143958868894602 0.106463878326996],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtNumOfFractions',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editFraction';

h22 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','30',...
'Style','edit',...
'Position',[0.219794344473008 0.228136882129278 0.161953727506427 0.0836501901140684],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editFraction_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editFraction_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editFraction',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtIsoCenter';

h23 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','IsoCenter in [mm]',...
'Style','text',...
'Position',[0.0282776349614396 0.330798479087452 0.201799485861183 0.091254752851711],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtIsoCenter',...
'UserData',[],...
'FontName','Helvetica',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editIsoCenter';

h24 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','0 0 0',...
'Style','edit',...
'Position',[0.219794344473008 0.338403041825095 0.161953727506427 0.0836501901140684],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editIsoCenter_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editIsoCenter_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editIsoCenter',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'checkIsoCenter';

h25 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Auto.',...
'Style','checkbox',...
'Value',1,...
'Position',[0.38560411311054 0.338403041825095 0.0809768637532133 0.091254752851711],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('checkIsoCenter_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','checkIsoCenter',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnRunSequencing';

h26 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Run Sequencing',...
'Style','radiobutton',...
'Position',[0.553984575835475 0.628020880324805 0.173521850899743 0.140684410646388],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnRunSequencing_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btnRunSequencing',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnRunDAO';

h27 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Run Direct Aperture Optimization',...
'Style','radiobutton',...
'Position',[0.553984575835475 0.32003608945028 0.294344473007712 0.140684410646388],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnRunDAO_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btnRunDAO',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtSequencing';

h28 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Stratification Levels',...
'Style','text',...
'Position',[0.553984575835475 0.502545595153702 0.179948586118252 0.102661596958175],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtSequencing',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'editSequencingLevel';

h29 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','7',...
'Style','edit',...
'Position',[0.58611825192802 0.449313655990204 0.0668380462724936 0.0836501901140685],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('editSequencingLevel_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('editSequencingLevel_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','editSequencingLevel',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'popUpMachine';

h30 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',['Generic';'       '],...
'Style','popupmenu',...
'Value',1,...
'Position',[0.219794344473008 0.418250950570342 0.161953727506427 0.114068441064639],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popUpMachine_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popUpMachine_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popUpMachine',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtMachine';

h31 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Machine',...
'Style','text',...
'Position',[0.0732647814910026 0.433460076045627 0.106683804627249 0.0950570342205323],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtMachine',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnSetTissue';

h32 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Set Tissue',...
'Position',[0.401028277634961 0.110266159695817 0.109254498714653 0.0874524714828897],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnSetTissue_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btnSetTissue',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popMenuBioOpt';

h33 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'none'; 'const_RBExD'; 'LEMIV_effect'; 'LEMIV_RBExD' },...
'Style','popupmenu',...
'Value',1,...
'Position',[0.219794344473008 0.0760456273764259 0.165809768637532 0.11787072243346],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popMenuBioOpt_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popMenuBioOpt_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popMenuBioOpt',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'text38';

h34 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Type of optimization',...
'Style','text',...
'Position',[0.0102827763496144 0.0988593155893536 0.201799485861183 0.091254752851711],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Interruptible','off',...
'Tag','text38',...
'UserData',[],...
'FontName','Helvetica',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobutton3Dconf';

h35 = uicontrol(...
'Parent',h12,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','3D conformal',...
'Style','radiobutton',...
'Position',[0.553224553224553 0.757869249394673 0.212121212121212 0.0847457627118644],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobutton3Dconf_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','radiobutton3Dconf',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel2';

h36 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'ShadowColor',get(0,'defaultuipanelShadowColor'),...
'Title',{  'Visualization'; '             ' },...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel2',...
'UserData',[],...
'Clipping','off',...
'Position',[0.00451321727917473 0.0460947503201024 0.430689877498388 0.203585147247119],...
'FontName','Helvetica',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupPlane';

h37 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'coronal'; 'sagital'; 'axial' },...
'Style','popupmenu',...
'Value',3,...
'Position',[0.465315808689303 0.582191780821918 0.113636363636364 0.143835616438356],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupPlane_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupPlane_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupPlane',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'sliderSlice';

h38 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'Slider' },...
'Style','slider',...
'Position',[0.134961439588689 0.796610169491525 0.167095115681234 0.096045197740113],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('sliderSlice_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'BusyAction','cancel',...
'Interruptible','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('sliderSlice_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','sliderSlice');

appdata = [];
appdata.lastValidTag = 'txtPlanSelection';

h39 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Plane Selection',...
'Style','text',...
'Position',[0.34258853596203 0.589041095890411 0.11969696969697 0.116438356164384],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtPlanSelection',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text9';

h40 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Slice Selection',...
'Style','text',...
'Position',[0.00909090909090909 0.808219178082192 0.113636363636364 0.116438356164384],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text9',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnContour';

h41 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot contour',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.780445969125214 0.733212341197822 0.169811320754717 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnContour_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnContour',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnDose';

h42 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot dose',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.780445969125214 0.466969147005444 0.169811320754717 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnDose_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnDose',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnIsoDoseLines';

h43 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot isolines',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.780445969125214 0.600907441016334 0.150943396226415 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnIsoDoseLines_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnIsoDoseLines',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtTypeOfPlot';

h44 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Type of plot',...
'Style','text',...
'Position',[0.343053173241852 0.793103448275862 0.108061749571183 0.124137931034483],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtTypeOfPlot',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupTypeOfPlot';

h45 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'intensity'; 'profile' },...
'Style','popupmenu',...
'Value',1,...
'Position',[0.465315808689303 0.801369863013699 0.113636363636364 0.143835616438356],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupTypeOfPlot_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupTypeOfPlot_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupTypeOfPlot',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtDisplayOption';

h46 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Display option',...
'Style','text',...
'Position',[0.34258853596203 0.39041095890411 0.128787878787879 0.102739726027397],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtDisplayOption',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupDisplayOption';

h47 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Please select ... ',...
'Style','popupmenu',...
'Value',1,...
'Position',[0.465315808689303 0.36986301369863 0.196969696969697 0.136986301369863],...
'BackgroundColor',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupDisplayOption_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupDisplayOption_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupDisplayOption',...
'FontWeight','bold');

appdata = [];
appdata.lastValidTag = 'txtBeamSelection';

h48 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Beam Selection',...
'Style','text',...
'Position',[0.00857632933104631 0.503448275862069 0.118353344768439 0.186206896551724],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtBeamSelection',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'sliderBeamSelection';

h49 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','SliderBeamSelection',...
'Style','slider',...
'Position',[0.134961439588689 0.542372881355932 0.167095115681234 0.096045197740113],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('sliderBeamSelection_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('sliderBeamSelection_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','sliderBeamSelection');

appdata = [];
appdata.lastValidTag = 'btnProfileType';

h50 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','lateral',...
'Position',[0.658025928757913 0.794520547945205 0.0863636363636364 0.157534246575343],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnProfileType_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btnProfileType',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text16';

h51 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','GoTo',...
'Style','text',...
'Position',[0.604795511376553 0.801369863013698 0.0454545454545454 0.123287671232877],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text16',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnDVH';

h52 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Show DVH/QI',...
'Position',[0.51413881748072 0.0677966101694915 0.123393316195373 0.129943502824859],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnDVH_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btnDVH',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnIsoDoseLinesLabels';

h53 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot isolines labels',...
'Style','radiobutton',...
'Position',[0.780445969125214 0.343557168784029 0.202401372212693 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnIsoDoseLinesLabels_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnIsoDoseLinesLabels',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'textOffset';

h54 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Offset',...
'Style','text',...
'Position',[0.00909090909090909 0.287104393008975 0.118181818181818 0.123287671232877],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','textOffset',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'sliderOffset';

h55 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','SliderOffset',...
'Style','slider',...
'Position',[0.134961439588689 0.271186440677966 0.167095115681234 0.096045197740113],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('sliderOffset_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('sliderOffset_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','sliderOffset');

appdata = [];
appdata.lastValidTag = 'radioBtnIsoCenter';

h56 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot iso center',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.780445969125214 0.205989110707804 0.169811320754717 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radioBtnIsoCenter_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radioBtnIsoCenter',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btn3Dview';

h57 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Open 3D-View',...
'Position',[0.595848595848596 0.578947368421053 0.148962148962149 0.157894736842105],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btn3Dview_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Enable','off',...
'Tag','btn3Dview',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnCT';

h58 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','plot CT',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.780445969125214 0.864791288566243 0.169811320754717 0.117241379310345],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnCT_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnCT',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'radiobtnPlan';

h59 = uicontrol(...
'Parent',h36,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','visualize plan / beams',...
'Style','radiobutton',...
'Value',1,...
'Position',[0.78021978021978 0.0736842105263158 0.2002442002442 0.115789473684211],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('radiobtnPlan_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','radiobtnPlan',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uitoolbar1';

h60 = uitoolbar(...
'Parent',h1,...
'Tag','uitoolbar1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Standard.FileOpen';
appdata.CallbackInUse = struct(...
    'ClickedCallback', 'matRadGUI(''toolbarLoad_ClickedCallback'',gcbo,[],guidata(gcbo))');
appdata.lastValidTag = 'toolbarLoad';

h61 = uipushtool(...
'Parent',h60,...
'Children',[],...
'BusyAction','cancel',...
'Interruptible','off',...
'Tag','toolbarLoad',...
'CData',mat{2},...
'ClickedCallback',@(hObject,eventdata)matRadGUI_export('toolbarLoad_ClickedCallback',hObject,eventdata,guidata(hObject)),...
'Separator','on',...
'TooltipString','Open File',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Standard.SaveFigure';
appdata.CallbackInUse = struct(...
    'ClickedCallback', 'matRadGUI(''toolbarSave_ClickedCallback'',gcbo,[],guidata(gcbo))');
appdata.lastValidTag = 'toolbarSave';

h62 = uipushtool(...
'Parent',h60,...
'Children',[],...
'BusyAction','cancel',...
'Interruptible','off',...
'Tag','toolbarSave',...
'CData',mat{3},...
'ClickedCallback',@(hObject,eventdata)matRadGUI_export('toolbarSave_ClickedCallback',hObject,eventdata,guidata(hObject)),...
'Separator','on',...
'TooltipString','Save Figure',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = [];
appdata.lastValidTag = 'uipushtool_screenshot';

h63 = uipushtool(...
'Parent',h60,...
'Children',[],...
'Tag','uipushtool_screenshot',...
'CData',mat{4},...
'ClickedCallback',@(hObject,eventdata)matRadGUI_export('uipushtool_screenshot_ClickedCallback',hObject,eventdata,guidata(hObject)),...
'TooltipString','Take a screenshot of the current dose or profile plot',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.ZoomIn';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'toolbarZoomIn';

h64 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','toolbarZoomIn',...
'CData',mat{5},...
'ClickedCallback','%default',...
'Separator','on',...
'TooltipString','Zoom In',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.ZoomOut';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'toolbarZoomOut';

h65 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','toolbarZoomOut',...
'CData',mat{6},...
'ClickedCallback','%default',...
'Separator','on',...
'TooltipString','Zoom Out',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.Pan';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'toolbarPan';

h66 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','toolbarPan',...
'CData',mat{7},...
'ClickedCallback','%default',...
'Separator','on',...
'TooltipString','Pan',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Exploration.DataCursor';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'toolbarCursor';

h67 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','toolbarCursor',...
'CData',mat{8},...
'ClickedCallback','%default',...
'Separator','on',...
'TooltipString','Data Cursor',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Annotation.InsertLegend';
appdata.CallbackInUse = struct(...
    'ClickedCallback', '%default');
appdata.lastValidTag = 'toolbarLegend';

h68 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','toolbarLegend',...
'CData',mat{9},...
'ClickedCallback','%default',...
'Separator','on',...
'TooltipString','Insert Legend',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.toolid = 'Annotation.InsertColorbar';
appdata.CallbackInUse = struct(...
    'ClickedCallback', 'matRadGUI(''uitoggletool8_ClickedCallback'',gcbo,[],guidata(gcbo))');
appdata.lastValidTag = 'uitoggletool8';

h69 = uitoggletool(...
'Parent',h60,...
'Children',[],...
'Tag','uitoggletool8',...
'CData',mat{10},...
'ClickedCallback',@(hObject,eventdata)matRadGUI_export('uitoggletool8_ClickedCallback',hObject,eventdata,guidata(hObject)),...
'Separator','on',...
'TooltipString','Insert Colorbar',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel3';

h70 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'ShadowColor',get(0,'defaultuipanelShadowColor'),...
'Title',{  'Objectives & constraints'; '                        '; '                        ' },...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel3',...
'UserData',[],...
'Clipping','off',...
'Position',[0.00451321727917473 0.257362355953905 0.430689877498388 0.259923175416133],...
'FontName','Helvetica',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel4';

h71 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'ShadowColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Title','Workflow',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel4',...
'UserData',[],...
'Clipping','off',...
'Position',[0.00451321727917473 0.810499359795134 0.430045132172792 0.170294494238156],...
'FontName','Helvetica',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text13';

h72 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Status:',...
'Style','text',...
'Position',[0.318250377073907 0.107438016528926 0.120663650075415 0.12396694214876],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text13',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtInfo';

h73 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','no data loaded',...
'Style','text',...
'Position',[0.414781297134238 0.0247933884297521 0.371040723981901 0.214876033057851],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtInfo',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnLoadMat';

h74 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Load  *.mat data',...
'Position',[0.151866151866152 0.810126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnLoadMat_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnLoadMat',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnCalcDose';

h75 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Calc. influence Mx',...
'Position',[0.35006435006435 0.810126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnCalcDose_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnCalcDose',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnOptimize';

h76 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Optimize',...
'Position',[0.544401544401541 0.810126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnOptimize_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnOptimize',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnLoadDicom';

h77 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Load DICOM',...
'Position',[0.151866151866152 0.60126582278481 0.177606177606178 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnLoadDicom_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnLoadDicom',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnRefresh';

h78 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Refresh',...
'Position',[0.0154440154440154 0.810126582278481 0.0849420849420849 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnRefresh_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnRefresh',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton_recalc';

h79 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Recalc',...
'Position',[0.543114543114543 0.60126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('pushbutton_recalc_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','pushbutton_recalc',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnSaveToGUI';

h80 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Save to GUI',...
'Position',[0.738738738738737 0.810126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnSaveToGUI_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnSaveToGUI',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btn_export';

h81 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Export',...
'Position',[0.74002574002574 0.60126582278481 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btn_export_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btn_export',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'importDoseButton';

h82 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Import Dose',...
'Position',[0.738738738738738 0.392405063291139 0.178893178893179 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('importDoseButton_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','importDoseButton',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton_importFromBinary';

h83 = uicontrol(...
'Parent',h71,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Import from Binary',...
'Position',[0.151866151866152 0.392405063291139 0.177606177606178 0.145569620253165],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('pushbutton_importFromBinary_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'TooltipString','Imports a patient data set from binary datafiles describing CT and segmentations',...
'Tag','pushbutton_importFromBinary',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text24';

h84 = uicontrol(...
'Parent',h1,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','min value:',...
'Style','text',...
'Position',[0.899701069855255 0.87145643693108 0.0420862177470107 0.0253576072821847],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text24',...
'FontSize',10,...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnSetIsoDoseLevels';

h85 = uicontrol(...
'Parent',h1,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Set IsoDose Levels',...
'Position',[0.910792951541851 0.814995131450828 0.071035242290749 0.0223953261927945],...
'BackgroundColor',[0.8 0.8 0.8],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnSetIsoDoseLevels_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnSetIsoDoseLevels',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel10';

h86 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'Title','Structure Visibilty',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel10',...
'Clipping','off',...
'Position',[0.896397105097545 0.175812743823147 0.0991189427312775 0.304291287386216],...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'legendTable';

h87 = uicontrol(...
'Parent',h86,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'Style','listbox',...
'Value',1,...
'Position',[0.02 0.01 0.97 0.98],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('legendTable_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('legendTable_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','legendTable');

appdata = [];
appdata.lastValidTag = 'uipanel11';

h88 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'Title','Viewing',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel11',...
'Clipping','off',...
'Position',[0.437782076079949 0.0460947503201025 0.451321727917473 0.842509603072983],...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'axesFig';

h89 = axes(...
'Parent',h88,...
'FontUnits',get(0,'defaultaxesFontUnits'),...
'Units',get(0,'defaultaxesUnits'),...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'CameraTarget',[0.5 0.5 0.5],...
'CameraTargetMode',get(0,'defaultaxesCameraTargetMode'),...
'CameraViewAngle',6.60861036031192,...
'CameraViewAngleMode',get(0,'defaultaxesCameraViewAngleMode'),...
'PlotBoxAspectRatio',[1 0.994838709677419 0.994838709677419],...
'PlotBoxAspectRatioMode',get(0,'defaultaxesPlotBoxAspectRatioMode'),...
'FontName','CMU Serif',...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'ColormapMode',get(0,'defaultaxesColormapMode'),...
'Alphamap',[0 0.0159 0.0317 0.0476 0.0635 0.0794 0.0952 0.1111 0.127 0.1429 0.1587 0.1746 0.1905 0.2063 0.2222 0.2381 0.254 0.2698 0.2857 0.3016 0.3175 0.3333 0.3492 0.3651 0.381 0.3968 0.4127 0.4286 0.4444 0.4603 0.4762 0.4921 0.5079 0.5238 0.5397 0.5556 0.5714 0.5873 0.6032 0.619 0.6349 0.6508 0.6667 0.6825 0.6984 0.7143 0.7302 0.746 0.7619 0.7778 0.7937 0.8095 0.8254 0.8413 0.8571 0.873 0.8889 0.9048 0.9206 0.9365 0.9524 0.9683 0.9841 1],...
'AlphamapMode',get(0,'defaultaxesAlphamapMode'),...
'GridLineStyle',get(0,'defaultaxesGridLineStyle'),...
'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
'XTickMode',get(0,'defaultaxesXTickMode'),...
'XTickLabel',{  '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1' },...
'XTickLabelMode',get(0,'defaultaxesXTickLabelMode'),...
'YTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1],...
'YTickMode',get(0,'defaultaxesYTickMode'),...
'YTickLabel',{  '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1' },...
'YTickLabelMode',get(0,'defaultaxesYTickLabelMode'),...
'Color',get(0,'defaultaxesColor'),...
'CameraMode',get(0,'defaultaxesCameraMode'),...
'DataSpaceMode',get(0,'defaultaxesDataSpaceMode'),...
'ColorSpaceMode',get(0,'defaultaxesColorSpaceMode'),...
'DecorationContainerMode',get(0,'defaultaxesDecorationContainerMode'),...
'ChildContainerMode',get(0,'defaultaxesChildContainerMode'),...
'XRulerMode',get(0,'defaultaxesXRulerMode'),...
'YRulerMode',get(0,'defaultaxesYRulerMode'),...
'ZRulerMode',get(0,'defaultaxesZRulerMode'),...
'AmbientLightSourceMode',get(0,'defaultaxesAmbientLightSourceMode'),...
'ButtonDownFcn',@(hObject,eventdata)matRadGUI_export('axesFig_ButtonDownFcn',hObject,eventdata,guidata(hObject)),...
'Tag','axesFig',...
'UserData',[],...
'Position',[0.0718390804597701 0.0354391371340524 0.902298850574712 0.929121725731895],...
'ActivePositionProperty','position',...
'ActivePositionPropertyMode',get(0,'defaultaxesActivePositionPropertyMode'),...
'LooseInset',[0.199881557553276 0.118695547917832 0.146067292058163 0.0809287826712493],...
'FontSize',9.63,...
'SortMethod','childorder',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h90 = get(h89,'title');

set(h90,...
'Parent',h89,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0 0 0],...
'ColorMode','auto',...
'Position',[0.500000554361651 1.00343482490272 0.5],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','Helvetica',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',2,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','Axes Title',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h91 = get(h89,'xlabel');

set(h91,...
'Parent',h89,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0.500000476837158 -0.0366861225689741 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','top',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h92 = get(h89,'ylabel');

set(h92,...
'Parent',h89,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[-0.0429806458688552 0.500000476837158 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',90,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10.593,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','bottom',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','back',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','on',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

h93 = get(h89,'zlabel');

set(h93,...
'Parent',h89,...
'Units','data',...
'FontUnits','points',...
'DecorationContainer',[],...
'DecorationContainerMode','auto',...
'Color',[0.15 0.15 0.15],...
'ColorMode','auto',...
'Position',[0 0 0],...
'PositionMode','auto',...
'String',blanks(0),...
'Interpreter','tex',...
'Rotation',0,...
'RotationMode','auto',...
'FontName','CMU Serif',...
'FontSize',10,...
'FontSizeMode','auto',...
'FontAngle','normal',...
'FontWeight','normal',...
'HorizontalAlignment','left',...
'HorizontalAlignmentMode','auto',...
'VerticalAlignment','middle',...
'VerticalAlignmentMode','auto',...
'EdgeColor','none',...
'LineStyle','-',...
'LineWidth',0.5,...
'BackgroundColor','none',...
'Margin',3,...
'Clipping','off',...
'Layer','middle',...
'LayerMode','auto',...
'FontSmoothing','on',...
'FontSmoothingMode','auto',...
'DisplayName',blanks(0),...
'IncludeRenderer','on',...
'IsContainer','off',...
'IsContainerMode','auto',...
'HelpTopicKey',blanks(0),...
'ButtonDownFcn',blanks(0),...
'BusyAction','queue',...
'Interruptible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} ,...
'DeleteFcn',blanks(0),...
'Tag',blanks(0),...
'HitTest','on',...
'PickableParts','visible',...
'PickablePartsMode','auto',...
'DimensionNames',{  'X' 'Y' 'Z' },...
'DimensionNamesMode','auto',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'Description','AxisRulerBase Label',...
'DescriptionMode','auto',...
'Visible','off',...
'Serializable','on',...
'HandleVisibility','off',...
'TransformForPrintFcnImplicitInvoke','on',...
'TransformForPrintFcnImplicitInvokeMode','auto');

appdata = [];
appdata.lastValidTag = 'uipanel12';

h94 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'Title','Info',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel12',...
'Clipping','off',...
'Position',[0.896276240708709 0.0448143405889885 0.0991189427312775 0.12932138284251],...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btnAbout';

h95 = uicontrol(...
'Parent',h94,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','About',...
'Position',[0.238095238095238 0.134831460674157 0.563492063492063 0.280898876404494],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Callback',@(hObject,eventdata)matRadGUI_export('btnAbout_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','btnAbout',...
'FontSize',7,...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text15';

h96 = uicontrol(...
'Parent',h94,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','v3.0.0',...
'Style','text',...
'Position',[0.227106227106227 0.752808988764045 0.523809523809524 0.191011235955056],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text15',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text31';

h97 = uicontrol(...
'Parent',h94,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','github.com/e0404/matRad',...
'Style','text',...
'Position',[0.0384615384615385 0.528089887640449 0.942307692307693 0.168539325842697],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text31',...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'uipanel_colormapOptions';

h98 = uipanel(...
'Parent',h1,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units',get(0,'defaultuipanelUnits'),...
'Title','Viewer Options',...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Tag','uipanel_colormapOptions',...
'Clipping','off',...
'Position',[0.896397105097545 0.484330299089727 0.0991189427312775 0.318660598179457],...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text_windowCenter';

h99 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Window Center:',...
'Style','text',...
'Position',[0.0466666666666666 0.682461750109027 0.673333333333333 0.0699999999999998],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','text_windowCenter',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'textDoseOpacity';

h100 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Dose opacity:',...
'Style','text',...
'Position',[0.0466666666666667 0.0706370831711431 0.847328244274809 0.0714285714285714],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','textDoseOpacity',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu_chooseColorData';

h101 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'None'; 'CT (ED)'; 'Dose' },...
'Style','popupmenu',...
'Value',1,...
'Position',[0.0486486486486487 0.899328859060403 0.940540540540541 0.11744966442953],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupmenu_chooseColorData_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupmenu_chooseColorData_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupmenu_chooseColorData');

appdata = [];
appdata.lastValidTag = 'slider_windowCenter';

h102 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'SliderStep',[0.01 0.05],...
'String','slider',...
'Style','slider',...
'Value',0.5,...
'Position',[0.0432432432432432 0.63758389261745 0.697297297297297 0.0536912751677853],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('slider_windowCenter_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('slider_windowCenter_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','slider_windowCenter');

appdata = [];
appdata.lastValidTag = 'text_windowWidth';

h103 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Window Width:',...
'Style','text',...
'Position',[0.0466666666666667 0.545761302394105 0.673333333333333 0.0700000000000001],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','text_windowWidth',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu_chooseColormap';

h104 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Choose Colormap...',...
'Style','popupmenu',...
'Value',1,...
'Position',[0.0362903225806452 0.278843516266481 0.939516129032258 0.0844686648501362],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupmenu_chooseColormap_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupmenu_chooseColormap_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupmenu_chooseColormap');

appdata = [];
appdata.lastValidTag = 'text_windowRange';

h105 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Range:',...
'Style','text',...
'Position',[0.0403225806451613 0.387807911050966 0.274193548387097 0.0708446866485015],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','text_windowRange',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'edit_windowRange';

h106 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','0 1',...
'Style','edit',...
'Position',[0.323863636363636 0.399846781328902 0.653409090909091 0.0707395498392283],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('edit_windowRange_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('edit_windowRange_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','edit_windowRange');

appdata = [];
appdata.lastValidTag = 'edit_windowCenter';

h107 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','0.5',...
'Style','edit',...
'Value',1,...
'Position',[0.767567567567568 0.63758389261745 0.205405405405405 0.0704697986577181],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('edit_windowCenter_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('edit_windowCenter_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','edit_windowCenter');

appdata = [];
appdata.lastValidTag = 'edit_windowWidth';

h108 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','1.0',...
'Style','edit',...
'Position',[0.772727272727273 0.518256759964609 0.204545454545455 0.0707395498392284],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('edit_windowWidth_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('edit_windowWidth_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','edit_windowWidth');

appdata = [];
appdata.lastValidTag = 'sliderOpacity';

h109 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'SliderStep',[0.01 0.05],...
'String','slider',...
'Style','slider',...
'Value',0.6,...
'Position',[0.147727272727273 0.0257234726688103 0.75 0.0546623794212219],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('sliderOpacity_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('sliderOpacity_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','sliderOpacity');

appdata = [];
appdata.lastValidTag = 'txtDoseOpacity0Indicator';

h110 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','0',...
'Style','text',...
'Position',[0.0466666666666666 0.00599285798906697 0.0810810810810811 0.072463768115942],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','txtDoseOpacity0Indicator',...
'UserData',[],...
'FontName','Helvetica',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtDoseOpacity1Indicator';

h111 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','right',...
'String','1',...
'Style','text',...
'Position',[0.8963482566536 0.00864864051690258 0.0810810810810811 0.072463768115942],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','txtDoseOpacity1Indicator',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text_windowPreset';

h112 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'HorizontalAlignment','left',...
'String','Window Presets N/A',...
'Style','text',...
'Position',[0.0540540540540541 0.842281879194631 0.697297297297297 0.0704697986577181],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'Children',[],...
'Tag','text_windowPreset',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu_windowPreset';

h113 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String',{  'Custom'; 'Full'; 'Abd/Med'; 'Head'; 'Liver'; 'Lung'; 'Spine'; 'Vrt/Bone' },...
'Style','popupmenu',...
'Value',1,...
'Position',[0.0486486486486487 0.73489932885906 0.940540540540541 0.11744966442953],...
'BackgroundColor',[1 1 1],...
'Callback',@(hObject,eventdata)matRadGUI_export('popupmenu_windowPreset_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Visible','off',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('popupmenu_windowPreset_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','popupmenu_windowPreset');

appdata = [];
appdata.lastValidTag = 'slider_windowWidth';

h114 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'SliderStep',[0.01 0.05],...
'String','slider',...
'Style','slider',...
'Value',1,...
'Position',[0.0454545454545455 0.51140507995425 0.698863636363636 0.0546623794212219],...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@(hObject,eventdata)matRadGUI_export('slider_windowWidth_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)matRadGUI_export('slider_windowWidth_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'Tag','slider_windowWidth');

appdata = [];
appdata.lastValidTag = 'checkbox_lockColormap';

h115 = uicontrol(...
'Parent',h98,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','Lock Settings',...
'Style','checkbox',...
'Position',[0.0486486486486487 0.151006711409396 0.940540540540541 0.0838926174496644],...
'BackgroundColor',[0.502 0.502 0.502],...
'Callback',@(hObject,eventdata)matRadGUI_export('checkbox_lockColormap_Callback',hObject,eventdata,guidata(hObject)),...
'Children',[],...
'Tag','checkbox_lockColormap',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text39';

h116 = uicontrol(...
'Parent',h1,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','max  value:',...
'Style','text',...
'Position',[0.901903713027061 0.85370611183355 0.0420862177470107 0.0245123537061118],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','text39',...
'FontSize',10,...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtMinVal';

h117 = uicontrol(...
'Parent',h1,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','-',...
'Style','text',...
'Position',[0.955789804908748 0.879908972691808 0.0271397105097546 0.0160598179453836],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtMinVal',...
'FontSize',10,...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'txtMaxVal';

h118 = uicontrol(...
'Parent',h1,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','normalized',...
'String','-',...
'Style','text',...
'Position',[0.955789804908748 0.863003901170351 0.0271397105097546 0.0177503250975293],...
'BackgroundColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
'Children',[],...
'Tag','txtMaxVal',...
'FontSize',10,...
'FontWeight','bold',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error(message('MATLAB:guide:StateFieldNotFound', gui_StateFields{ i }, gui_Mfile));
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % MATRADGUI_EXPORT
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % MATRADGUI_EXPORT(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallback(gui_State, varargin{:})
    % MATRADGUI_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % MATRADGUI_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~ishghandle(fig,'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || (isscalar(fig)&&isprop(fig,'GUIDEFigure'));
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishghandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI MATLAB code file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    % Call version of openfig that accepts 'auto' option"
    gui_hFigure = matlab.hg.internal.openfigLegacy(name, singleton, visible);  
%     %workaround for CreateFcn not called to create ActiveX
%         peers=findobj(findall(allchild(gui_hFigure)),'type','uicontrol','style','text');    
%         for i=1:length(peers)
%             if isappdata(peers(i),'Control')
%                 actxproxy(peers(i));
%             end            
%         end
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallback(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishghandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end


