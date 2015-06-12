function varargout = matRad_importDicomGUI(varargin)
% MATRAD_IMPORTDICOMGUI MATLAB code for matRad_importDicomGUI.fig
%      MATRAD_IMPORTDICOMGUI, by itself, creates a new MATRAD_IMPORTDICOMGUI or raises the existing
%      singleton*.
%
%      H = MATRAD_IMPORTDICOMGUI returns the handle to a new MATRAD_IMPORTDICOMGUI or the handle to
%      the existing singleton*.
%
%      MATRAD_IMPORTDICOMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRAD_IMPORTDICOMGUI.M with the given input arguments.
%
%      MATRAD_IMPORTDICOMGUI('Property','Value',...) creates a new MATRAD_IMPORTDICOMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matRad_importDicomGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matRad_importDicomGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matRad_importDicomGUI

% Last Modified by GUIDE v2.5 11-Jun-2015 18:41:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @matRad_importDicomGUI_OpeningFcn, ...
    'gui_OutputFcn',  @matRad_importDicomGUI_OutputFcn, ...
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


% --- Executes just before matRad_importDicomGUI is made visible.
function matRad_importDicomGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRad_importDicomGUI (see VARARGIN)

% Choose default command line output for matRad_importDicomGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes matRad_importDicomGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matRad_importDicomGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse_button.
function patDir = browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiwait(warndlg('Choose the input directory'));
patDir = uigetdir('', 'Choose the input directory...');
if patDir ~= 0
    patDir = [patDir '\'];
    handles.dir_path_field.String = patDir;
    % Update handles structure
    guidata(hObject, handles);
    %     [fileList, patient_listbox] = matRad_scanDicomImportFolder_h(handles.dir_path_field.String);
    %     if iscell(patient_listbox)
    %         handles.fileList =  fileList;
    %         handles.patient_listbox.String = patient_listbox;
    %         % Update handles structure
    %         guidata(hObject, handles);
    %     end
    scan(hObject, eventdata, handles)
end

function scan(hObject, eventdata, handles)
[fileList, patient_listbox] = matRad_scanDicomImportFolder(handles.dir_path_field.String);
if iscell(patient_listbox)
    handles.fileList =  fileList;
    handles.patient_listbox.String = patient_listbox;
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on selection change in patient_listbox.
function patient_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to patient_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns patient_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from patient_listbox

if ~isempty(hObject.String)
    % enable Import button
    handles.import_button.Enable = 'on';
    
    % handles.filelist:
    %   1. Filepath
    %   2. Modality
    %   3. PatientID
    %   4. SeriesUID
    %   5. SeriesNumber
    %   9. res_x
    %   10. res_y
    %   11. res_z
    selected_patient = handles.patient_listbox.String(handles.patient_listbox.Value);
    if handles.SeriesUID_radiobutton.Value == 1
        % this gets a list of ct series for this patient
        handles.ct_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4));
        % this gets a list of rtss series for this patient
        handles.rtss_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4));
        % this gets a resolution for this patient
        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), handles.ct_listbox.String),9));
        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), handles.ct_listbox.String),10));
        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), handles.ct_listbox.String),11));
    else
        handles.ct_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5));
        handles.rtss_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),5));
        res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), handles.ct_listbox.String),9));
        res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), handles.ct_listbox.String),10));
        res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), handles.ct_listbox.String),11));
    end
    handles.resx_edit.String = res_x;
    handles.resy_edit.String = res_y;
    handles.resz_edit.String = res_z;
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function patient_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to patient_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ct_listbox.
function ct_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to ct_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ct_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ct_listbox


% --- Executes during object creation, after setting all properties.
function ct_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ct_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rtss_listbox.
function rtss_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to rtss_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rtss_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rtss_listbox


% --- Executes during object creation, after setting all properties.
function rtss_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rtss_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in import_button.
function import_button_Callback(hObject, eventdata, handles)
% hObject    handle to import_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selected_patient = handles.patient_listbox.String(handles.patient_listbox.Value);
selected_ctseries = handles.ct_listbox.String(handles.ct_listbox.Value);
selected_rtseries = handles.rtss_listbox.String(handles.rtss_listbox.Value);

if handles.SeriesUID_radiobutton.Value == 1
    files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,4), selected_ctseries),:);
    files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,4), selected_rtseries),:);
else
    files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,5), selected_ctseries),:);
    files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,5), selected_rtseries),:);
end

files.resx = str2double(handles.resx_edit.String);
files.resy = str2double(handles.resy_edit.String);
files.resz = str2double(handles.resz_edit.String);
matRad_importDicom(files);


% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);


% --- Executes on button press in rescan_button.
function rescan_button_Callback(hObject, eventdata, handles)
% hObject    handle to rescan_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function dir_path_field_Callback(hObject, eventdata, handles)
% hObject    handle to dir_path_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dir_path_field as text
%        str2double(get(hObject,'String')) returns contents of dir_path_field as a double

patDir = handles.dir_path_field.String;
if patDir(end) ~= '\';
    patDir = [patDir '\'];
    handles.dir_path_field.String = patDir;
    % Update handles structure
    guidata(hObject, handles);
end
scan(hObject, eventdata, handles);


% --- Executes on button press in SeriesUID_radiobutton.
function SeriesUID_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to SeriesUID_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value == 1
    handles.SeriesNumber_radiobutton.Value = 0;
else
    hObject.Value = 1;
    handles.SeriesNumber_radiobutton.Value = 0;
end
selected_patient = handles.patient_listbox.String(handles.patient_listbox.Value);
handles.ct_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4));
handles.rtss_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4));
guidata(hObject, handles);

% --- Executes on button press in SeriesNumber_radiobutton.
function SeriesNumber_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to SeriesNumber_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1
    handles.SeriesUID_radiobutton.Value = 0;
else
    hObject.Value = 1;
    handles.SeriesUID_radiobutton.Value = 0;
end
selected_patient = handles.patient_listbox.String(handles.patient_listbox.Value);
handles.ct_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5));
handles.rtss_listbox.String = unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),5));
guidata(hObject, handles);



function resx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to resx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resx_edit as text
%        str2double(get(hObject,'String')) returns contents of resx_edit as a double


% --- Executes during object creation, after setting all properties.
function resx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resy_edit_Callback(hObject, eventdata, handles)
% hObject    handle to resy_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resy_edit as text
%        str2double(get(hObject,'String')) returns contents of resy_edit as a double


% --- Executes during object creation, after setting all properties.
function resy_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resy_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function resz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to resz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resz_edit as text
%        str2double(get(hObject,'String')) returns contents of resz_edit as a double


% --- Executes during object creation, after setting all properties.
function resz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
