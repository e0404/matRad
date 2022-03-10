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

% Last Modified by GUIDE v2.5 02-Jun-2017 00:45:04

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
    'gui_OpeningFcn', @matRad_importDicomGUI_OpeningFcn, ...
    'gui_OutputFcn',  @matRad_importDicomGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

% if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
% else
%     gui_mainfcn(gui_State, varargin{:});
% end
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

axes(handles.axesMatRadLogo)
[im, ~, alpha] = imread('matrad_logo.png');
q = image(im);
axis equal off
set(q, 'AlphaData', alpha);
% show dkfz logo
axes(handles.axesDKFZLogo)
[im, ~, alpha] = imread('DKFZ_Logo.png');
p = image(im);
axis equal off
set(p, 'AlphaData', alpha);
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
    patDir = [patDir filesep];
    %handles.dir_path_field.String = patDir;
    set(handles.dir_path_field,'String',patDir);
    % Update handles structure
    guidata(hObject, handles);
    scan(hObject, eventdata, handles)
end

function scan(hObject, eventdata, handles)
[fileList, patient_listbox] = matRad_scanDicomImportFolder(get(handles.dir_path_field,'String'));
if iscell(patient_listbox)
    handles.fileList =  fileList;
    %handles.patient_listbox.String = patient_listbox;
    set(handles.patient_listbox,'String',patient_listbox,'Value',1);
    guidata(hObject, handles);
end

% --- Executes on selection change in patient_listbox.
function patient_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to patient_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns patient_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from patient_listbox

if ~isempty(get(hObject,'String'))
    % enable Import button
    set(handles.import_button,'Enable','on');
    
    % handles.filelist:
    %   1. Filepath
    %   2. Modality
    %   3. PatientID
    %   4. SeriesUID
    %   5. SeriesNumber
    %   9. res_x
    %   10. res_y
    %   11. res_z
    %   12. detailed dose description - currently not in use for GUI user
    patient_listbox = get(handles.patient_listbox,'String');
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
    % this gets a list of rtss series for this patient
    set(handles.rtseries_listbox,'Value',1); % set dummy value to one
    set(handles.rtseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4));
    % this gets a list of rt plan series for this patient
    set(handles.rtplan_listbox,'Value',[]); % set dummy value to none
    set(handles.rtplan_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTPLAN') & strcmp(handles.fileList(:,3), selected_patient),4));
    % this gets a list of dose series for this patient
    set(handles.doseseries_listbox,'Value',[]); % set dummy value to none
    set(handles.doseseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTDOSE') & strcmp(handles.fileList(:,3), selected_patient),4));
    % selectedDose

    if get(handles.SeriesUID_radiobutton,'Value') == 1
        % this gets a list of ct series for this patient
        set(handles.ctseries_listbox,'Value',1); % set dummy value to one
        set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4)));
        
        selectedDoseSeriesString = get(handles.doseseries_listbox,'String');
        % this gets a resolution for this patient
        selectedCtSeriesString = get(handles.ctseries_listbox,'String');
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    else
        set(handles.ctseries_listbox,'Value',1); % set dummy value to one
        set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5)));
        selectedCtSeriesString = get(handles.ctseries_listbox,'String');
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    end
    set(handles.resx_edit,'String',res_x);
    set(handles.resy_edit,'String',res_y);
    if numel(res_z) > 1
        set(handles.resz_edit,'String','not equi');
    else
        set(handles.resz_edit,'String',res_z);
    end
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


% --- Executes on selection change in ctseries_listbox.
function ctseries_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to ctseries_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ctseries_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ctseries_listbox


% --- Executes during object creation, after setting all properties.
function ctseries_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ctseries_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in rtseries_listbox.
function rtseries_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to rtseries_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rtseries_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rtseries_listbox


% --- Executes during object creation, after setting all properties.
function rtseries_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rtseries_listbox (see GCBO)
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

patient_listbox = get(handles.patient_listbox,'String');
ctseries_listbox = get(handles.ctseries_listbox,'String');
rtplan_listbox = get(handles.rtplan_listbox,'String');
doseseries_listbox = get(handles.rtplan_listbox,'String');
if ~isempty(patient_listbox)
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
end
if ~isempty(ctseries_listbox)
    selected_ctseries = ctseries_listbox(get(handles.ctseries_listbox,'Value'));
end
if ~isempty(rtplan_listbox)
    selected_rtplan = rtplan_listbox(get(handles.rtplan_listbox,'Value'));
end

if get(handles.SeriesUID_radiobutton,'Value') == 1
    files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,4), selected_ctseries),:);
    
    %files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
    %    strcmp(handles.fileList(:,4), selected_rtseries),:);
else
      files.ct = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
        strcmp(handles.fileList(:,5), selected_ctseries) & strcmp(handles.fileList(:,2),'CT'),:);    
    
    %files.rtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & ...
    %    strcmp(handles.fileList(:,5), selected_rtseries),:);
end

allRtss = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTSTRUCT'),:);
if ~isempty(allRtss)
    files.rtss = allRtss(get(handles.rtseries_listbox,'Value'),:);
else
    files.rtss = [];
end

files.resx = str2double(get(handles.resx_edit,'String'));
files.resy = str2double(get(handles.resy_edit,'String'));
% check if valid assignment is made when z slices are not equi-distant
if strcmp(get(handles.resz_edit,'String'),'not equi')
    msgbox('Ct data not sliced equi-distantly in z direction! Chose uniform resolution.', 'Error','error');
    return;
else
    files.resz = str2double(get(handles.resz_edit,'String'));
end
% selected RT Plan
rtplan = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTPLAN'),:);
if ~isempty(rtplan) && ~isempty(get(handles.rtplan_listbox,'Value'))
    files.rtplan = rtplan(get(handles.rtplan_listbox,'Value'),:);
end

% selected RT Dose
rtdose = handles.fileList(strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,2),'RTDOSE'),:);
if ~isempty(rtdose) && ~isempty(get(handles.doseseries_listbox,'Value'))
    selectedRtDose   = get(handles.doseseries_listbox,'String');
    selectedRtDoseIx = NaN*ones(1,numel(selectedRtDose));
    for i = 1:numel(selectedRtDose)
        selectedRtDoseIx(i) = find(strcmp(rtdose(:,4),selectedRtDose{i}));
    end
    files.rtdose = rtdose(selectedRtDoseIx,:);
end

% check if we should use the dose grid resolution
files.useDoseGrid = get(handles.checkbox3,'Value');

% dicomMetaBool: store complete DICOM information and patientName or not
dicomMetaBool = logical(get(handles.checkPatientName,'Value'));
[ct,cst,pln,stf,resultGUI] = matRad_importDicom(files, dicomMetaBool);


% Save file
matRadFileName = [files.ct{1,3} '.mat']; % use default from dicom
[FileName,PathName] = uiputfile('*','Save as...',matRadFileName);
if ischar(FileName)
    %Hack to get non-empty variables
    varNames = {'ct','cst','pln','stf','resultGUI'};
    ix = cellfun(@(s) ~isempty(evalin('caller',s)),varNames);
    varNames = varNames(ix);

    FileName = fullfile(PathName,FileName);

    % save all non-empty variables imported
    save(FileName,'-v7.3',varNames{:});
end


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

patDir = get(handles.dir_path_field,'String');
if patDir(end) ~= filesep;
    patDir = [patDir filesep];
    set(handles.dir_path_field,'String',patDir);
    guidata(hObject, handles);
end
scan(hObject, eventdata, handles);


% --- Executes on button press in SeriesUID_radiobutton.
function SeriesUID_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to SeriesUID_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') == 1
    set(handles.SeriesNumber_radiobutton,'Value',0);
else
    set(hObject,'Value',1);
    set(handles.SeriesNumber_radiobutton,'Value',0);
end
if isfield(handles, 'fileList')
    patient_listbox = get(handles.patient_listbox,'String');
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
    set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),4)));
    set(handles.rtseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTSTRUCT') & strcmp(handles.fileList(:,3), selected_patient),4)));
    set(handles.doseseries_listbox,'String',handles.fileList(strcmp(handles.fileList(:,2), 'RTDOSE') & strcmp(handles.fileList(:,3), selected_patient),4));
    set(handles.rtplan_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'RTPLAN') & strcmp(handles.fileList(:,3), selected_patient),4)));
else
    fprintf('No patient loaded, so just switching default display option to SeriesUID. \n');
end
guidata(hObject, handles);

% --- Executes on button press in SeriesNumber_radiobutton.
function SeriesNumber_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to SeriesNumber_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    set(handles.SeriesUID_radiobutton,'Value',0);
else
    set(hObject,'Value',1);
    set(handles.SeriesUID_radiobutton,'Value',0);
end
if isfield(handles, 'fileList')
    patient_listbox = get(handles.patient_listbox,'String');
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
    set(handles.ctseries_listbox,'String',unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient),5)));
else
    fprintf('No patient loaded, so just switching default display option to SeriesNumber. \n');
end
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

% 
% % --- Executes on selection change in ctseries_listbox.
% function ctseries_listbox_Callback(hObject, eventdata, handles)
% % hObject    handle to ctseries_listbox (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns ctseries_listbox contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from ctseries_listbox
% 
% 
% % --- Executes during object creation, after setting all properties.
% function ctseries_listbox_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to ctseries_listbox (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% % --- Executes on selection change in rtseries_listbox.
% function rtseries_listbox_Callback(hObject, eventdata, handles)
% % hObject    handle to rtseries_listbox (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns rtseries_listbox contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from rtseries_listbox
% 
% 
% % --- Executes during object creation, after setting all properties.
% function rtseries_listbox_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to rtseries_listbox (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on selection change in doseseries_listbox.
function doseseries_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to doseseries_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns doseseries_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from doseseries_listbox

if ~isempty(get(hObject,'Value'))
    set(handles.checkbox3,'Enable','on');
else
    set(handles.checkbox3,'Value',0);
    set(handles.checkbox3,'Enable','off');
    % retrieve and display resolution for DICOM ct cube
    patient_listbox = get(handles.patient_listbox,'String');
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
    selectedCtSeriesString = get(handles.ctseries_listbox,'String');
    if get(handles.SeriesUID_radiobutton,'Value') == 1
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    else
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    end
    set(handles.resx_edit,'String',res_x);
    set(handles.resy_edit,'String',res_y);
    if numel(res_z) > 1
        set(handles.resz_edit,'String','not equi');
    else
        set(handles.resz_edit,'String',res_z);
    end
    
end

% --- Executes during object creation, after setting all properties.
function doseseries_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to doseseries_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in rtplan_listbox.
function rtplan_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to rtplan_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
if ~isempty(get(hObject,'Value')) && numel(get(hObject,'Value')) == 1

    selectedPlan = contents{get(hObject,'Value')};
    % point at plan in listbox
    selectedPlanLoc = strcmp(handles.fileList(:,4),selectedPlan);
    
    % show only the doses corresponding to the plan
    corrDoses = [handles.fileList{selectedPlanLoc,13}];
    numOfDoses = numel(corrDoses);
    corrDosesLoc = zeros(size(handles.fileList(:,1),1),1);
    for j = 1:numOfDoses
        if ~isnan(corrDoses{j})
            corrDosesLoc = corrDosesLoc | strcmp(handles.fileList(:,4),corrDoses{j});
        end
    end
    
    if sum(corrDosesLoc) == 0
        warndlg('no rt dose file directly associated to plan file. showing all rt dose files.');
        corrDosesLoc = strcmp(handles.fileList(:,2),'RTDOSE');
    end

    set(handles.doseseries_listbox,'Value',[]); % set dummy value to one
    set(handles.doseseries_listbox,'String',handles.fileList(corrDosesLoc,4));
    
    % disable checkbox for use dose grid is currently checked
    if get(handles.checkbox3,'Value') == 1
        set(handles.checkbox3,'Value',0);
        checkbox3_Callback(handles.checkbox3,[], handles);
    end
    set(handles.checkbox3,'Enable','off');

    
elseif numel(get(hObject,'Value')) >=2
    warning('More than one RTPLAN selected. Unsetting selection ...');
    patient_listbox_Callback(hObject, eventdata, handles)
else
    patient_listbox_Callback(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function rtplan_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rtplan_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkPatientName.
function checkPatientName_Callback(hObject, eventdata, handles)
% hObject    handle to checkPatientName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%A = get(hObject,'Value');

% Hint: get(hObject,'Value') returns toggle state of checkPatientName
%guidata(hObject, handles);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

if get(hObject,'Value')
    set(handles.resx_edit,'Enable', 'off');
    set(handles.resy_edit,'Enable', 'off');
    set(handles.resz_edit,'Enable', 'off');
    % retrieve and display resolution for DICOM dose cube
    doseFilesInList = get(handles.doseseries_listbox,'String');
    selectedDoseFiles = get(handles.doseseries_listbox,'Value');
    if isempty(selectedDoseFiles)
        set(hObject,'Value',0)
        errordlg('no dose file selected');
        return;
    end
    for i = 1:numel(selectedDoseFiles)
        selectedDoseFile = doseFilesInList{selectedDoseFiles(i)};
        if verLessThan('matlab','9')
            dicomDoseInfo = dicominfo(handles.fileList{find(strcmp(handles.fileList(:,4),selectedDoseFile)),1});
        else
            dicomDoseInfo = dicominfo(handles.fileList{find(strcmp(handles.fileList(:,4),selectedDoseFile)),1},'UseDictionaryVR',true);
        end
        res_x{i} = dicomDoseInfo.PixelSpacing(1);
        res_y{i} = dicomDoseInfo.PixelSpacing(2);
        res_z{i} = dicomDoseInfo.SliceThickness;
    end

    if numel(unique(cell2mat(res_x)))*numel(unique(cell2mat(res_y)))*numel(unique(cell2mat(res_z))) ~= 1
        set(handles.checkbox3,'Value',0);
        warndlg('Different resolutions in dose file(s)');
        set(handles.resx_edit,'Enable', 'on');
        set(handles.resy_edit,'Enable', 'on');
        set(handles.resz_edit,'Enable', 'on');
    else
        set(handles.resx_edit,'String',num2str(res_x{1}));
        set(handles.resy_edit,'String',num2str(res_y{1}));
        set(handles.resz_edit,'String',num2str(res_z{1}));
    end
    
else
    set(handles.resx_edit,'Enable', 'on');
    set(handles.resy_edit,'Enable', 'on');
    set(handles.resz_edit,'Enable', 'on');
    % retrieve and display resolution for DICOM ct cube
    patient_listbox = get(handles.patient_listbox,'String');
    selected_patient = patient_listbox(get(handles.patient_listbox,'Value'));
    selectedCtSeriesString = get(handles.ctseries_listbox,'String');
    if get(handles.SeriesUID_radiobutton,'Value') == 1
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,4), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    else
        if ~isempty(selectedCtSeriesString)
            res_x = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),9));
            res_y = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),10));
            res_z = unique(handles.fileList(strcmp(handles.fileList(:,2), 'CT') & strcmp(handles.fileList(:,3), selected_patient) & strcmp(handles.fileList(:,5), selectedCtSeriesString{get(handles.ctseries_listbox,'Value')}),11));
        else
            res_x = NaN; res_y = NaN; res_z = NaN;
        end
    end
    set(handles.resx_edit,'String',res_x);
    set(handles.resy_edit,'String',res_y);
    if numel(res_z) > 1
        set(handles.resz_edit,'String','not equi');
    else
        set(handles.resz_edit,'String',res_z);
    end
    
end
