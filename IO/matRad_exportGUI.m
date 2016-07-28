function varargout = matRad_exportGUI(varargin)
% MATRAD_EXPORTGUI MATLAB code for matRad_exportGUI.fig
%      MATRAD_EXPORTGUI, by itself, creates a new MATRAD_EXPORTGUI or raises the existing
%      singleton*.
%
%      H = MATRAD_EXPORTGUI returns the handle to a new MATRAD_EXPORTGUI or the handle to
%      the existing singleton*.
%
%      MATRAD_EXPORTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATRAD_EXPORTGUI.M with the given input arguments.
%
%      MATRAD_EXPORTGUI('Property','Value',...) creates a new MATRAD_EXPORTGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matRad_exportGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matRad_exportGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help matRad_exportGUI

% Last Modified by GUIDE v2.5 07-Jul-2016 14:50:05

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
                   'gui_OpeningFcn', @matRad_exportGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @matRad_exportGUI_OutputFcn, ...
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

% --- Executes just before matRad_exportGUI is made visible.
function matRad_exportGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matRad_exportGUI (see VARARGIN)

% Choose default command line output for matRad_exportGUI
handles.output = hObject;

%Fills structure export table
if evalin('base','exist(''cst'',''var'')') == 1
    cst = evalin( 'base', 'cst' );
    tableData = cell(numel(cst(:,2)),2);
    tableData(:,2) = cst(:,2);
    tableData(:,1) = {true};
else
    tableData = cell(0);
    set(handles.checkbox_CT,'Enable','off');
end
set(handles.uitable_vois,'data',tableData);

%Fills result cubes export table
if evalin('base','exist(''resultGUI'',''var'')')
    result = evalin( 'base', 'resultGUI' );
    cubeNames = fieldnames(result);
    cubeIx = 1;
    for f = 1:numel(cubeNames)
        if ndims(result.(cubeNames{f})) < 3
            continue;
        end
        cubes{cubeIx} = cubeNames{f};
        cubeIx = cubeIx + 1;
    end    
    numCubes = cubeIx - 1;
    tableData = cell(numCubes,2);
    tableData(:,2) = cubes;
    tableData(:,1) = {true};
else
    tableData = cell(0);
    set(handles.checkbox_dose,'Enable','off');
end
set(handles.uitable_doseCubes,'data',tableData);


% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes matRad_exportGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = matRad_exportGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the btn_cancel flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to btn_cancel the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

%{
handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

set(handles.density, 'String', handles.metricdata.density);
set(handles.volume,  'String', handles.metricdata.volume);
set(handles.mass, 'String', 0);

set(handles.unitgroup, 'SelectedObject', handles.english);

set(handles.text4, 'String', 'lb/cu.in');
set(handles.text5, 'String', 'cu.in');
set(handles.text6, 'String', 'lb');

% Update handles structure
guidata(handles.figure1, handles);
%}

% --- Executes on button press in checkbox_CT.
function checkbox_CT_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveCT = get(hObject,'Value');

%Show the VOI-table only if we want to save a CT
if (saveCT)
    set(handles.uitable_vois,'Visible', 'on', 'Enable','on');
else
    set(handles.uitable_vois,'Visible', 'off', 'Enable','off');
end    


% --- Executes on selection change in listbox_vois.
function uitable_vois_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vois contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vois


% --- Executes during object creation, after setting all properties.
function uitable_vois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_dose.
function checkbox_dose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Show the Result-table only if we want to save dose cubes
saveDose = get(hObject,'Value');
if (saveDose)
    set(handles.uitable_doseCubes,'Visible', 'on', 'Enable','on');

else
    set(handles.uitable_doseCubes,'Visible', 'off', 'Enable','off');
    %set(handles.uitable_vois,'data',cell(0));
end    


% --- Executes on selection change in listbox_dose.
function uitable_doseCubes_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_dose contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_dose


% --- Executes during object creation, after setting all properties.
function listbox_dose_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_dose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_dir_export_browse.
function exportDir = pushbutton_dir_export_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dir_export_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

exportDir = uigetdir('', 'Choose the input directory...');
if exportDir ~= 0
    exportDir = [exportDir filesep];
    set(handles.edit_dir_export,'String',exportDir);
    % Update handles structure
    guidata(hObject, handles);
end


function exportDir = edit_dir_export_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dir_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dir_export as text
%        str2double(get(hObject,'String')) returns contents of edit_dir_export as a double

exportDir = get(handles.edit_dir_export,'String');

%Add filesperator
if exportDir(end) ~= filesep;
    exportDir = [exportDir filesep];
end

%Check if the user specified an existing directory
if ~exist(exportDir,'dir')
    warndlg(['Folder ' exportDir ' does not exist!']);
    exportDir = '';
end    
set(handles.edit_dir_export,'String',exportDir);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_dir_export_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dir_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    

    
end


% --- Executes on selection change in popupmenu_extension.
function popupmenu_extension_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_extension contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_extension


% --- Executes during object creation, after setting all properties.
function popupmenu_extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%These sets up the available extensions
extensions{1} = '*.nrrd';
extensions{2} = '*.vtk';
extensions{3} = '*.mha';
set(hObject,'String',extensions);

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_export.
function btn_export_Callback(hObject, eventdata, handles)
% hObject    handle to btn_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get the export dir
exportDir = get(handles.edit_dir_export,'String');

%Sanity check
if numel(exportDir) == 0
    errordlg('No Export folder selected!');
    return;
elseif ~exist(exportDir,'dir')
    errordlg(['Folder ' exportDir ' does not exist!']);
    return;
else
    %Add file separator if necessary
    if exportDir(end) ~= filesep;
        exportDir = [exportDir filesep];
    end 
end

%Get the file extension
extensionIndex = get(handles.popupmenu_extension,'Value');
extensions = get(handles.popupmenu_extension,'String');
extension = extensions{extensionIndex};
extension = extension(2:end);

saveCT = get(handles.checkbox_CT,'Value');
saveResults = get(handles.checkbox_dose,'Value');

%%Prepare for export
%If we export CT, try to create a subdirectory for VOIs
if (saveCT)
    voiDir = [exportDir '/vois/'];
    if ~exist(voiDir,'dir')
        if ~mkdir(voiDir)
            warndlg('Could not create subfolder for VOI masks. Masks will be stored in base folder.');
            voiDir = exportDir;
        end
    end
end
%If we export results, try to create a subdirectory for VOIs
if (saveResults)
    resultDir = [exportDir '/results/'];
    if ~exist(resultDir,'dir')
        if ~mkdir(resultDir)
            warndlg('Could not create subfolder for resulting dose cubes. Cubes will be stored in base folder.');
            resultDir = exportDir;
        end
    end
end

%prepare metadata
ct = evalin('base','ct');

metadata.resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];
metadata.compress = get(handles.checkbox_compress,'Value');

%Check if we have position information
if isfield(ct,'dicomInfo')
    if isfield(ct.dicomInfo,'ImagePositionPatient')
       metadata.imageOrigin = ct.dicomInfo.ImagePositionPatient;       
    end    
end

%This is only for the waitbar to get the number of cubes you wanna save
if (saveCT)
    numExportCubes = 1;
    voiNames = get(handles.uitable_vois,'Data');
    voiIndices = find([voiNames{:,1}] == true);
    numExportCubes = numExportCubes + numel(voiIndices);
else
    numExportCubes = 0;
end
if saveResults
   cubeNames = get(handles.uitable_doseCubes,'data');
   cubeIndices = find([cubeNames{:,1}] == true);
   numExportCubes = numExportCubes + numel(cubeIndices);
end

%Give an error if nothing was selected
if numExportCubes == 0
    errordlg('No data was selected for export!');
    return;
end       

currentCube = 0;

hWaitbar = waitbar(0,'Exporting...','WindowStyle', 'modal');
cleanUp = onCleanup(@() close(hWaitbar));

%CT and Mask export
if saveCT   
    %Export the CT (ED suffix to clarify it is not in HU)
    currentCube = currentCube + 1;
    waitbar(currentCube/numExportCubes,hWaitbar,['Exporting CT (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
    
    matRad_writeCube(fullfile(exportDir,['CT_ED' extension]),ct.cube{1},'double',metadata);
    
    %Export VOI masks
    cst = evalin('base','cst');
    
    for voiIx = voiIndices
        %Waitbar
        currentCube = currentCube + 1;
        waitbar(currentCube/numExportCubes,hWaitbar,['Exporting Segmentation Mask (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
        
        %Get the index list
        voiRow = find(strcmp(voiNames{voiIx,2},cst(:,2)));
        voiIndexList = cst{voiRow,4}{1};
        %Set up the full mask
        voiMask = zeros(ct.cubeDim);
        voiMask(voiIndexList) = 1;
        %Export...
        matRad_writeCube(fullfile(voiDir,[voiNames{voiIx,2} extension]),voiMask,'uint8',metadata);
    end
    
end

%Results Export
if saveResults
    results = evalin('base','resultGUI');    
    cubeNames = get(handles.uitable_doseCubes,'data');
       
    for cubeIx = cubeIndices
        %Export
        currentCube = currentCube + 1;
        waitbar(currentCube/numExportCubes,hWaitbar,['Exporting Results (' num2str(currentCube) '/' num2str(numExportCubes) ') ...']);
        matRad_writeCube(fullfile(resultDir,[cubeNames{cubeIx,2} extension]),results.(cubeNames{cubeIx,2}),'double',metadata);
    end
end

close(handles.figure1);
    

% --- Executes on button press in btn_cancel.
function btn_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

% --- Executes during object creation, after setting all properties.
function btn_cancel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to btn_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in checkbox_compress.
function checkbox_compress_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_compress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_compress
