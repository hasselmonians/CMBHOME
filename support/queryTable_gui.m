function varargout = queryTable_gui(varargin)
% QUERYTABLE_GUI MATLAB code for queryTable_gui.fig
%      QUERYTABLE_GUI, by itself, creates a new QUERYTABLE_GUI or raises the existing
%      singleton*.
%
%      H = QUERYTABLE_GUI returns the handle to a new QUERYTABLE_GUI or the handle to
%      the existing singleton*.
%
%      QUERYTABLE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUERYTABLE_GUI.M with the given input arguments.
%
%      QUERYTABLE_GUI('Property','Value',...) creates a new QUERYTABLE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before queryTable_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to queryTable_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help queryTable_gui

% Last Modified by GUIDE v2.5 19-Feb-2013 14:09:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @queryTable_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @queryTable_gui_OutputFcn, ...
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


% --- Executes just before queryTable_gui is made visible.
function queryTable_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to queryTable_gui (see VARARGIN)

% Choose default command line output for queryTable_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes queryTable_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = queryTable_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in queryButton.
function queryButton_Callback(hObject, eventdata, handles)
% hObject    handle to queryButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume




function rat_name_Callback(hObject, eventdata, handles)
% hObject    handle to rat_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_name as text
%        str2double(get(hObject,'String')) returns contents of rat_name as a double


% --- Executes during object creation, after setting all properties.
function rat_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rat_notes_Callback(hObject, eventdata, handles)
% hObject    handle to rat_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_notes as text
%        str2double(get(hObject,'String')) returns contents of rat_notes as a double


% --- Executes during object creation, after setting all properties.
function rat_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rat_implant1_Callback(hObject, eventdata, handles)
% hObject    handle to rat_implant1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_implant1 as text
%        str2double(get(hObject,'String')) returns contents of rat_implant1 as a double


% --- Executes during object creation, after setting all properties.
function rat_implant1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_implant1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rat_implant2_Callback(hObject, eventdata, handles)
% hObject    handle to rat_implant2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_implant2 as text
%        str2double(get(hObject,'String')) returns contents of rat_implant2 as a double


% --- Executes during object creation, after setting all properties.
function rat_implant2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_implant2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rat_depth1_Callback(hObject, eventdata, handles)
% hObject    handle to rat_depth1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_depth1 as text
%        str2double(get(hObject,'String')) returns contents of rat_depth1 as a double


% --- Executes during object creation, after setting all properties.
function rat_depth1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_depth1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rat_depth2_Callback(hObject, eventdata, handles)
% hObject    handle to rat_depth2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rat_depth2 as text
%        str2double(get(hObject,'String')) returns contents of rat_depth2 as a double


% --- Executes during object creation, after setting all properties.
function rat_depth2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rat_depth2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in grp.
function grp_Callback(hObject, eventdata, handles)
% hObject    handle to grp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns grp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from grp


% --- Executes during object creation, after setting all properties.
function grp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tetrode_notes_Callback(hObject, eventdata, handles)
% hObject    handle to tetrode_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tetrode_notes as text
%        str2double(get(hObject,'String')) returns contents of tetrode_notes as a double


% --- Executes during object creation, after setting all properties.
function tetrode_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tetrode_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tetrode_num2_Callback(hObject, eventdata, handles)
% hObject    handle to tetrode_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tetrode_num2 as text
%        str2double(get(hObject,'String')) returns contents of tetrode_num2 as a double


% --- Executes during object creation, after setting all properties.
function tetrode_num2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tetrode_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tetrode_num1_Callback(hObject, eventdata, handles)
% hObject    handle to tetrode_num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tetrode_num1 as text
%        str2double(get(hObject,'String')) returns contents of tetrode_num1 as a double


% --- Executes during object creation, after setting all properties.
function tetrode_num1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tetrode_num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to session_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_filepath as text
%        str2double(get(hObject,'String')) returns contents of session_filepath as a double


% --- Executes during object creation, after setting all properties.
function session_filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_notes_Callback(hObject, eventdata, handles)
% hObject    handle to session_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_notes as text
%        str2double(get(hObject,'String')) returns contents of session_notes as a double


% --- Executes during object creation, after setting all properties.
function session_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_date1_Callback(hObject, eventdata, handles)
% hObject    handle to session_date1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_date1 as text
%        str2double(get(hObject,'String')) returns contents of session_date1 as a double


% --- Executes during object creation, after setting all properties.
function session_date1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_date1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_date2_Callback(hObject, eventdata, handles)
% hObject    handle to session_date2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_date2 as text
%        str2double(get(hObject,'String')) returns contents of session_date2 as a double


% --- Executes during object creation, after setting all properties.
function session_date2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_date2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_experimenter_Callback(hObject, eventdata, handles)
% hObject    handle to session_experimenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_experimenter as text
%        str2double(get(hObject,'String')) returns contents of session_experimenter as a double


% --- Executes during object creation, after setting all properties.
function session_experimenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_experimenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_manipulation_Callback(hObject, eventdata, handles)
% hObject    handle to session_manipulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_manipulation as text
%        str2double(get(hObject,'String')) returns contents of session_manipulation as a double


% --- Executes during object creation, after setting all properties.
function session_manipulation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_manipulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_task_Callback(hObject, eventdata, handles)
% hObject    handle to session_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_task as text
%        str2double(get(hObject,'String')) returns contents of session_task as a double


% --- Executes during object creation, after setting all properties.
function session_task_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_task (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ts_notes_Callback(hObject, eventdata, handles)
% hObject    handle to ts_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_notes as text
%        str2double(get(hObject,'String')) returns contents of ts_notes as a double


% --- Executes during object creation, after setting all properties.
function ts_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ts_depth1_Callback(hObject, eventdata, handles)
% hObject    handle to ts_depth1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_depth1 as text
%        str2double(get(hObject,'String')) returns contents of ts_depth1 as a double


% --- Executes during object creation, after setting all properties.
function ts_depth1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_depth1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ts_depth2_Callback(hObject, eventdata, handles)
% hObject    handle to ts_depth2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_depth2 as text
%        str2double(get(hObject,'String')) returns contents of ts_depth2 as a double


% --- Executes during object creation, after setting all properties.
function ts_depth2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_depth2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ts_anatomy_Callback(hObject, eventdata, handles)
% hObject    handle to ts_anatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ts_anatomy as text
%        str2double(get(hObject,'String')) returns contents of ts_anatomy as a double


% --- Executes during object creation, after setting all properties.
function ts_anatomy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ts_anatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_notes_Callback(hObject, eventdata, handles)
% hObject    handle to cell_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_notes as text
%        str2double(get(hObject,'String')) returns contents of cell_notes as a double


% --- Executes during object creation, after setting all properties.
function cell_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_num2_Callback(hObject, eventdata, handles)
% hObject    handle to cell_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_num2 as text
%        str2double(get(hObject,'String')) returns contents of cell_num2 as a double


% --- Executes during object creation, after setting all properties.
function cell_num2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_num1_Callback(hObject, eventdata, handles)
% hObject    handle to cell_num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_num1 as text
%        str2double(get(hObject,'String')) returns contents of cell_num1 as a double


% --- Executes during object creation, after setting all properties.
function cell_num1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_morphology_Callback(hObject, eventdata, handles)
% hObject    handle to cell_morphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_morphology as text
%        str2double(get(hObject,'String')) returns contents of cell_morphology as a double


% --- Executes during object creation, after setting all properties.
function cell_morphology_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_morphology (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cell_type_Callback(hObject, eventdata, handles)
% hObject    handle to cell_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_type as text
%        str2double(get(hObject,'String')) returns contents of cell_type as a double


% --- Executes during object creation, after setting all properties.
function cell_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_taskNotes_Callback(hObject, eventdata, handles)
% hObject    handle to session_taskNotes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_taskNotes as text
%        str2double(get(hObject,'String')) returns contents of session_taskNotes as a double


% --- Executes during object creation, after setting all properties.
function session_taskNotes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_taskNotes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function session_enclosure_Callback(hObject, eventdata, handles)
% hObject    handle to session_enclosure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session_enclosure as text
%        str2double(get(hObject,'String')) returns contents of session_enclosure as a double


% --- Executes during object creation, after setting all properties.
function session_enclosure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session_enclosure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epoch_label_Callback(hObject, eventdata, handles)
% hObject    handle to epoch_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch_label as text
%        str2double(get(hObject,'String')) returns contents of epoch_label as a double


% --- Executes during object creation, after setting all properties.
function epoch_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epoch_notes_Callback(hObject, eventdata, handles)
% hObject    handle to epoch_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epoch_notes as text
%        str2double(get(hObject,'String')) returns contents of epoch_notes as a double


% --- Executes during object creation, after setting all properties.
function epoch_notes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epoch_notes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ce_rm.
function ce_rm_Callback(hObject, eventdata, handles)
% hObject    handle to ce_rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ce_rm


% --- Executes on button press in ce_acorr.
function ce_acorr_Callback(hObject, eventdata, handles)
% hObject    handle to ce_acorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ce_acorr


% --- Executes on button press in includeIDs.
function includeIDs_Callback(hObject, eventdata, handles)
% hObject    handle to includeIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of includeIDs
