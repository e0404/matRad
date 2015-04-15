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

% Last Modified by GUIDE v2.5 15-Apr-2015 20:04:26

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

% Update handles structure
guidata(hObject, handles);

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
