function varargout = matRad_importGUI(varargin)
% MATRAD_IMPORTGUI MATLAB code for matRad_importGUI.fig
%      MATRAD_IMPORTGUI, by itself, creates a new MATRAD_IMPORTGUI or raises the existing
%      singleton*.
%
%      H = MATRAD_IMPORTGUI returns the handle to a new MATRAD_IMPORTGUI or the handle to
%      the existing singleton*.
%
%      MATRAD_IMPORTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRAD_IMPORTGUI.M with the given input arguments.
%
%      MATRAD_IMPORTGUI('Property','Value',...) creates a new MATRAD_IMPORTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matRad_importGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matRad_importGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matRad_importGUI

% Last Modified by GUIDE v2.5 09-Aug-2018 15:18:30

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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matRad_importGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @matRad_importGUI_OutputFcn, ...
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


% --- Executes just before matRad_importGUI is made visible.
function matRad_importGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRad_importGUI (see VARARGIN)

% Choose default command line output for matRad_importGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes matRad_importGUI wait for user response (see UIRESUME)
% uiwait(handles.figure_importDialog);


% --- Outputs from this function are returned to the command line.
function varargout = matRad_importGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_ctPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ctPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ctPath as text
%        str2double(get(hObject,'String')) returns contents of edit_ctPath as a double


% --- Executes during object creation, after setting all properties.
function edit_ctPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ctPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ctPath.
function pushbutton_ctPath_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ctPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[importCTFile,importCTPath,~] = uigetfile({'*.nrrd', 'NRRD-Files'}, 'Choose the CT file...');
if importCTFile ~= 0
    set(handles.edit_ctPath,'String',fullfile(importCTPath,importCTFile));
    % Update handles structure
    guidata(hObject, handles);
end



function listbox_maskPaths_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_maskPaths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of listbox_maskPaths as text
%        str2double(get(hObject,'String')) returns contents of listbox_maskPaths as a double


% --- Executes during object creation, after setting all properties.
function listbox_maskPaths_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_maskPaths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'value',[],'max',2,'min',0,'String',cell(0));


% --- Executes on button press in pushbutton_addMaskPaths.
function pushbutton_addMaskPaths_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addMaskPaths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[importMaskFile,importMaskPath,~] = uigetfile({'*.nrrd', 'NRRD-Files'}, 'Choose the binary mask files...','MultiSelect','on');
if ~isempty(importMaskFile)
    if ~iscell(importMaskFile)
        tmpName = importMaskFile;
        importMaskFile = cell(1);
        importMaskFile{1} = tmpName;
    end
    importMaskFile = cellfun(@(filename) fullfile(importMaskPath,filename),importMaskFile,'UniformOutput',false);
    entries = get(handles.listbox_maskPaths,'String');
    newEntries = [entries importMaskFile];
    set(handles.listbox_maskPaths,'String',newEntries);
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_import.
function pushbutton_import_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ctFile = get(handles.edit_ctPath,'String');
maskFiles = get(handles.listbox_maskPaths,'String');

if isempty(ctFile) || isempty(maskFiles)
    errordlg('Please sepecify a CT and at least one mask!');
end

convertHU = get(handles.checkbox_huConvert,'Value');

if convertHU
    [ct,cst] = matRad_importPatient(ctFile,maskFiles,get(handles.edit_hlut,'String'));
else
    [ct,cst] = matRad_importPatient(ctFile,maskFiles);
end

cst = showCheckDialog(cst);

assignin('base', 'ct', ct);
assignin('base', 'cst', cst);

delete(handles.figure_importDialog);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure_importDialog);


% --- Executes on button press in pushbutton_addMaskFolders.
function pushbutton_addMaskFolders_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addMaskFolders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
importMaskPath = uigetdir('./', 'Choose the folder containing binary mask files...');
importMaskPath = [importMaskPath filesep];
if ~isempty(importMaskPath)
    entries = get(handles.listbox_maskPaths,'String');
    newEntries = [entries cellstr(importMaskPath)];
    set(handles.listbox_maskPaths,'String',newEntries);
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on key press with focus on listbox_maskPaths and none of its controls.
function listbox_maskPaths_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_maskPaths (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if isequal(eventdata.Key,'delete') || isequal(eventdata.Key,'backspace')
    selectIndex = get(hObject,'value');
    entries = get(hObject,'String');
    if numel(entries) == 0
        return;
    end
    entries(selectIndex) = [];
    if selectIndex > numel(entries) && selectIndex > 1
        selectIndex = selectIndex - 1;
    end
    set(hObject,'String',entries,'value',selectIndex);
end

% --- Creates a Dialog for the final adaptations to VOIs and CT conversion
function cst = showCheckDialog(cst)

handle = dialog('Position', [100 100 400 250],'WindowStyle','modal','Name','Confirm Segmentations');

%Create Table
hTable = uitable('Parent',handle,'Units','normal','Position',[0.1 0.2 0.8 0.8]);
hTable.Data = cst(:,2:3);
hTable.ColumnName = {'Name','Type'};
hTable.ColumnWidth = {150,'auto'};
hTable.RowName = [];
hTable.ColumnEditable = [true true];
hTable.ColumnFormat = {'char',{'TARGET', 'OAR', 'IGNORED'}};

%Create Button
hButton = uicontrol(handle,'Style','pushbutton','String','Confirm','Units','normal','Position',[0.7 0.05 0.2 0.1],'Callback','uiresume(gcbf)');%{@pushbutton_confirm_vois_callback});
try
    uiwait(handle);
    cst(:,2:3) = hTable.Data(:,:);
catch
    warning('Closed checkdialog without confirmation! Using default cst information!');
end
delete(handle);


% --- Executes on button press in checkbox_huConvert.
function checkbox_huConvert_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_huConvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_huConvert

checked = get(hObject,'Value');

if checked
    fieldState = 'on';
else
    fieldState = 'off';
end


set(handles.edit_hlut,'Enable',fieldState);
set(handles.pushbutton_hlutFile,'Enable',fieldState);


function edit_hlut_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hlut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hlut as text
%        str2double(get(hObject,'String')) returns contents of edit_hlut as a double


% --- Executes during object creation, after setting all properties.
function edit_hlut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hlut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_hlutFile.
function pushbutton_hlutFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_hlutFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[importHLUTFile,importHLUTPath,~] = uigetfile({'*.hlut', 'matRad HLUT-Files'}, 'Choose the HLUT file...');
if importHLUTFile ~= 0
    set(handles.edit_hlut,'String',fullfile(importHLUTPath,importHLUTFile));
    % Update handles structure
    guidata(hObject, handles);
end
