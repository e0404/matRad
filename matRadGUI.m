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
    addpath(fullfile(matRadRootDir,'tools'))
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

function handles = resetGUI(hObject, handles, varargin)
% enable opengl software rendering to visualize linewidths properly
if ispc
  opengl software
elseif ismac
  % opengl is not supported
end


if ~isdeployed
    currFolder = fileparts(mfilename('fullpath'));
    addpath(fullfile(currFolder,'plotting'));
    addpath(fullfile(currFolder,['plotting' filesep 'colormaps']));
else
    currFolder = [];
end

% Choose default command line output for matRadGUI
handles.output = hObject;
%show matrad logo
axes(handles.axesLogo)
[im, ~, alpha] = imread([currFolder filesep 'dicomImport' filesep 'matrad_logo.png']);
f = image(im);
axis equal off
set(f, 'AlphaData', alpha);
% show dkfz logo
axes(handles.axesDKFZ)
[im, ~, alpha] = imread([currFolder filesep 'dicomImport' filesep 'DKFZ_Logo.png']);
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
        matRadRootDir = fileparts(mfilename('fullpath'));
        if ~isdeployed
            addpath(fullfile(matRadRootDir,'dicomImport'));
        end
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
     if ismember('pln',AllVarNames)  && handles.State > 0
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

% set slice slider
if handles.State > 0
    if evalin('base','exist(''pln'',''var'')')
        currPln = evalin('base','pln');
        if handles.plane == 1
            currSlice = ceil(currPln.propStf.isoCenter(1,2)/ct.resolution.x);
        elseif handles.plane == 2
            currSlice = ceil(currPln.propStf.isoCenter(1,1)/ct.resolution.y);
        elseif handles.plane == 3
            currSlice = ceil(currPln.propStf.isoCenter(1,3)/ct.resolution.z);
        end 
    else % no pln -> no isocenter -> use middle
        currSlice = ceil(ct.cubeDim(handles.plane)/2);
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

% --- Executes just before matRadGUI is made visible.
function matRadGUI_OpeningFcn(hObject, ~, handles, varargin) 
%#ok<*DEFNU> 
%#ok<*AGROW>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRadGUI (see VARARGIN)

% variable to check whether GUI is opened or just refreshed / new data
% loaded, since resetGUI needs to distinguish at one point

handles.initialGuiStart = true;

%If devMode is true, error dialogs will include the full stack trace of the error
%If false, only the basic error message is shown (works for errors that
%handle the MException object)
handles.devMode = true;

handles = resetGUI(hObject, handles);

%% parse variables from base workspace
AllVarNames = evalin('base','who');
handles.AllVarNames = AllVarNames;
try
    if  ismember('ct',AllVarNames) &&  ismember('cst',AllVarNames)
        ct  = evalin('base','ct');
        cst = evalin('base','cst');
        cst = setCstTable(handles,cst);
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
%guidata(findobj('Name','matRadGUI'), handles);
UpdatePlot(handles);

Update
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
    
catch
    handles = showWarning(handles,'LoadMatFileFnc: Could not load *.mat file');
    guidata(hObject,handles);
    UpdateState(handles);
    UpdatePlot(handles);
    return
end

try
    cst = setCstTable(handles,cst);
    handles.TableChanged = false;
    set(handles.popupTypeOfPlot,'Value',1);
    cst = matRad_computeVoiContoursWrapper(cst,ct);

    assignin('base','ct',ct);
    assignin('base','cst',cst);
    handles.State = 1;
catch
    handles = showError(handles,'LoadMatFileFnc: Could not load selected data');
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
    if ~isdeployed
        matRadRootDir = fileparts(mfilename('fullpath'));
        addpath(fullfile(matRadRootDir,'dicomImport'))
    end
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
    if ~isdeployed
        matRadRootDir = fileparts(mfilename('fullpath'));
        addpath(fullfile(matRadRootDir,'IO'))
    end
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
        set(handles.vmcFlag,'Value',0);
        set(handles.vmcFlag,'Enable','on')

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
        set(handles.vmcFlag,'Value',0);
        set(handles.vmcFlag,'Enable','off')
        
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
        set(handles.vmcFlag,'Value',0);
        set(handles.vmcFlag,'Enable','off')        
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
    
    pause(0.1);
    uiTable_CellEditCallback(hObject,[],handles);
    pause(0.3);

    %% get cst from table
    if ~getCstTable(handles);
        return
    end
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
        if get(handles.vmcFlag,'Value') == 0
            dij = matRad_calcPhotonDose(evalin('base','ct'),stf,pln,evalin('base','cst'));
        elseif get(handles.vmcFlag,'Value') == 1
            if ~isdeployed
                dij = matRad_calcPhotonDoseVmc(evalin('base','ct'),stf,pln,evalin('base','cst'));
            else
                error('VMC++ not available in matRad standalone application');
            end
        end
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
            handles.SelectedDisplayOption = 'physicalDose';
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
        handles.dispWindow{ctIx,2} = [min(ct.cubeHU{:}(:)) max(ct.cubeHU{:}(:))];
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
    set(handles.radiobtnPlan,'enable','on');
    axesFig3D = handles.axesFig3D;
    fig3D = handles.fig3D;    
else
    set(handles.radiobtnPlan,'value',0,'enable','off');
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

if ~isempty(pln)
    set(axesFig3D,'XTick',0:50:1000);
    set(axesFig3D,'YTick',0:50:1000);
    set(axesFig3D,'ZTick',0:50:1000);
    set(axesFig3D,'XTickLabel',0:50:1000);
    set(axesFig3D,'YTickLabel',0:50:1000);
    set(axesFig3D,'ZTickLabel',0:50:1000);
    xlabel(axesFig3D,'x [mm]','FontSize',defaultFontSize)
    ylabel(axesFig3D,'y [mm]','FontSize',defaultFontSize)
    zlabel(axesFig3D,'z [mm]','FontSize',defaultFontSize)
    %title(['axial plane z = ' num2str(ct.resolution.z*slice) ' [mm]'],'FontSize',defaultFontSize)
    title(axesFig3D,'3D view');
else
    xlabel(axesFig3D,'x [voxels]','FontSize',defaultFontSize)
    ylabel(axesFig3D,'y [voxels]','FontSize',defaultFontSize)
    zlabel(axesFig3D,'z [voxels]','FontSize',defaultFontSize)
    %title('axial plane','FontSize',defaultFontSize)
    title(axesFig3D,'matRad 3D view');
end

%if get(handles.radioBtnIsoCenter,'Value') == 1 && get(handles.popupTypeOfPlot,'Value') == 1 && ~isempty(pln)
%    hIsoCenterCross = matRad_plotIsoCenterMarker(handles.axesFig,pln,ct,plane,slice);
%end

% the following line ensures the plotting order (optional)
% set(gca,'Children',[AxesHandlesCT_Dose hIsoCenterCross AxesHandlesIsoDose  AxesHandlesVOI ]);

%set axis ratio
ratios = [1 1 1]; %[1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
ratios = ratios([2 1 3]);
set(axesFig3D,'DataAspectRatioMode','manual');
set(axesFig3D,'DataAspectRatio',ratios./max(ratios));

set(axesFig3D,'Ydir','reverse');

upperLimits = double(ct.cubeDim).*[ct.resolution.y ct.resolution.x ct.resolution.z];
set(axesFig3D,'xlim',[1 upperLimits(1)],'ylim',[1 upperLimits(2)],'zlim',[1 upperLimits(3)]);

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
    btnTableSave_Callback([],[],handles);


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
        [resultGUIcurrentRun,ipoptInfo] = matRad_fluenceOptimization(matRad_collapseDij(evalin('base','dij')),evalin('base','cst'),pln);
        resultGUIcurrentRun.w = resultGUIcurrentRun.w * ones(evalin('base','dij.totalNumOfBixels'),1);
        resultGUIcurrentRun.wUnsequenced = resultGUIcurrentRun.w;
    else
        if pln.propOpt.runDAO
        if ~matRad_checkForConnectedBixelRows(evalin('base','stf'))
            error('disconnetced dose influence data in BEV - run dose calculation again with consistent settings');
        end
        end
        
        [resultGUIcurrentRun,ipoptInfo] = matRad_fluenceOptimization(evalin('base','dij'),evalin('base','cst'),pln);
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
        CheckIpoptStatus(ipoptInfo,'Fluence')
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
       [resultGUI,ipoptInfo] = matRad_directApertureOptimization(evalin('base','dij'),evalin('base','cst'),...
           resultGUI.apertureInfo,resultGUI,pln);
       assignin('base','resultGUI',resultGUI);
       % check IPOPT status and return message for GUI user
       CheckIpoptStatus(ipoptInfo,'DAO');      
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
handles.dispWindow{3,1} = []; handles.dispWindow{3,2} = [];
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

% displays the cst in the GUI
function cst = setCstTable(handles,cst)

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

columnname = {'VOI name','VOI type','priority','obj. / const.','penalty','dose', 'EUD','volume','robustness'};

AllObjectiveFunction = {'square underdosing','square overdosing','square deviation', 'mean', 'EUD',...
                        'min dose constraint','max dose constraint',...
                        'min mean dose constraint','max mean dose constraint',...
                        'min EUD constraint','max EUD constraint',...
                        'max DVH constraint','min DVH constraint',...
                        'max DVH objective' ,'min DVH objective'};

columnformat = {cst(:,2)',{'OAR','TARGET'},'numeric',...
       AllObjectiveFunction,...
       'numeric','numeric','numeric','numeric',{'none','WC','prob'}};
   
numOfObjectives = 0;
for i = 1:size(cst,1)
    if ~isempty(cst{i,6})
        numOfObjectives = numOfObjectives + numel(cst{i,6});
    end
end

dimArr = [numOfObjectives size(columnname,2)];
data = cell(dimArr);
data(:,6) = {''};
Counter = 1;
for i = 1:size(cst,1)
   
   if strcmp(cst(i,3),'IGNORED')~=1
       for j=1:size(cst{i,6},1)
       %VOI
       data{Counter,1}  = cst{i,2};
       %VOI Type
       data{Counter,2}  = cst{i,3};
       %Priority
       data{Counter,3}  = cst{i,5}.Priority;
       %Objective Function
       objFunc          = cst{i,6}(j).type;
       data{Counter,4}  = objFunc;
       % set Penalty
       data{Counter,5}  = cst{i,6}(j).penalty;
       data{Counter,6}  = cst{i,6}(j).dose;
       data{Counter,7}  = cst{i,6}(j).EUD;
       data{Counter,8}  = cst{i,6}(j).volume;
       data{Counter,9}  = cst{i,6}(j).robustness;
       
       Counter = Counter +1;
       end
   end
   
end

set(handles.uiTable,'ColumnName',columnname);
set(handles.uiTable,'ColumnFormat',columnformat);
set(handles.uiTable,'ColumnEditable',[true true true true true true true true true true]);
set(handles.uiTable,'Data',data);


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
getCstTable(handles);
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
function CheckIpoptStatus(info,OptCase) 
      
if info.status == 0
    statusmsg = 'solved';  
elseif info.status == 1
    statusmsg = 'solved to acceptable level';          
elseif info.status == 2
    statusmsg = 'infeasible problem detected';           
elseif info.status == 3    
    statusmsg = 'search direction too small';             
elseif info.status == 4 
    statusmsg = 'diverging iterates';     
elseif info.status == 5
    statusmsg = 'user requested stop';     
elseif info.status == -1        
    statusmsg = 'maximum number of iterations';     
elseif info.status == -2    
    statusmsg = 'restoration phase failed';     
elseif info.status == -3         
    statusmsg = 'error in step computation';     
elseif info.status == -4
    statusmsg = 'maximum CPU time exceeded';     
elseif info.status == -10        
    statusmsg = 'not enough degrees of freedom';     
elseif info.status == -11    
    statusmsg = 'invalid problem definition';     
elseif info.status == -12   
    statusmsg = 'invalid option';     
elseif info.status == -13        
    statusmsg = 'invalid number detected';     
elseif info.status == -100
    statusmsg = 'unrecoverable exception';     
elseif info.status == -101        
    statusmsg = 'non-IPOPT exception thrown';     
elseif info.status == -102    
    statusmsg = 'insufficient memory';     
elseif info.status == -199
    statusmsg = 'IPOPT internal error';     
else
    statusmsg = 'IPOPT returned no status';
end
    
if info.status == 0 || info.status == 1
    status = 'none';
else
    status = 'warn';
end

msgbox(['IPOPT finished with status ' num2str(info.status) ' (' statusmsg ')'],'IPOPT',status,'modal');

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
        cst = setCstTable(handles,cst);
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
        if get(handles.vmcFlag,'Value') == 0
            dij = matRad_calcPhotonDose(ct,stf,pln,cst);
        elseif get(handles.vmcFlag,'Value') == 1
            dij = matRad_calcPhotonDoseVmc(evalin('base','ct'),stf,pln,evalin('base','cst'));
        end
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

% --- Executes on button press in vmcFlag.
function vmcFlag_Callback(hObject, eventdata, handles)
% hObject    handle to vmcFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vmcFlag


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
    tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
elseif ~handles.VOIPlotFlag(idx)
    handles.VOIPlotFlag(idx) = true;
    tmpString{idx} = ['<html><table border=0 ><TR><TD bgcolor=',clr,' width="18"><center>&#10004;</center></TD><TD>',cst{idx,2},'</TD></TR> </table></html>'];
end
set(handles.legendTable,'String',tmpString);

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
    matRadRootDir = fileparts(mfilename('fullpath'));
    if ~isdeployed
        addpath(fullfile(matRadRootDir,'IO'))
    end
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
    if ~isdeployed
        matRadRootDir = fileparts(mfilename('fullpath'));
        addpath(fullfile(matRadRootDir,'IO'))
    end
    
    %call the gui
    uiwait(matRad_importGUI);
    
    %Check if we have the variables in the workspace
    if evalin('base','exist(''cst'',''var'')') == 1 && evalin('base','exist(''ct'',''var'')') == 1
        cst = evalin('base','cst');
        ct = evalin('base','ct');
        setCstTable(handles,cst);
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

[filename, pathname] = uiputfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'; '*.fig','MATLAB figure file'},'Save current view','./screenshot.png');

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

%% Callbacks & Functions for color setting
function UpdateColormapOptions(handles)

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

target = get(event_obj,'Target');

%Get GUI data (maybe there is another way?)
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
