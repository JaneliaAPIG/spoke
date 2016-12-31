function varargout = SpikeGrid(varargin)
% SPIKEGRID MATLAB code for SpikeGrid.fig
%      SPIKEGRID, by itself, creates a new SPIKEGRID or raises the existing
%      singleton*.
%
%      H = SPIKEGRID returns the handle to a new SPIKEGRID or the handle to
%      the existing singleton*.
%
%      SPIKEGRID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKEGRID.M with the given input arguments.
%
%      SPIKEGRID('Property','Value',...) creates a new SPIKEGRID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpikeGrid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpikeGrid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpikeGrid

% Last Modified by GUIDE v2.5 28-Dec-2016 20:52:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpikeGrid_OpeningFcn, ...
                   'gui_OutputFcn',  @SpikeGrid_OutputFcn, ...
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


% --- Executes just before SpikeGrid is made visible.
function SpikeGrid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpikeGrid (see VARARGIN)

% Choose default command line output for SpikeGrid
handles.output = hObject;

% Update handles structure
handles.pcFreqAdjustProps = most.gui.control.PropertyTable(handles.tblFreqAdjustProps);
handles.pcFreqAdjustProps.alphabetize = false;
set(handles.tblFreqAdjustProps,'UserData',handles.pcFreqAdjustProps); % xxx this does not get set automatically by Controller probably b/c there is no "propbinding"

handles.pcRareAdjustProps = most.gui.control.PropertyTable(handles.tblRareAdjustProps);
handles.pcRareAdjustProps.alphabetize = false;
set(handles.tblRareAdjustProps,'UserData',handles.pcRareAdjustProps); % xxx this does not get set automatically by Controller probably b/c there is no "propbinding"

most.gui.AdvancedPanelToggler.init(hObject,handles.tbShowAllProps,16);


guidata(hObject, handles);

% UIWAIT makes SpikeGrid wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpikeGrid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbStartOrStop.
function pbStartOrStop_Callback(hObject, eventdata, handles)
% hObject    handle to pbStartOrStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pbStartOrStop
handles.hController.changeRunning();

% --- Executes on button press in cbGlobalMeanSubtract.
function cbGlobalMeanSubtract_Callback(hObject, eventdata, handles)
% hObject    handle to cbGlobalMeanSubtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbGlobalMeanSubtract
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes on button press in cbFilterData.
function cbFilterData_Callback(hObject, eventdata, handles)
% hObject    handle to cbFilterData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbFilterData
handles.hController.filterEnable = logical(get(hObject,'Value'));


function etThresholdVal_Callback(hObject, eventdata, handles)
% hObject    handle to etThresholdVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etThresholdVal as text
%        str2double(get(hObject,'String')) returns contents of etThresholdVal as a double
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function etThresholdVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etThresholdVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbDecThreshold.
function pbDecThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pbDecThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hController.stepThresholdVal(-1);


% --- Executes on button press in pbIncThreshold.
function pbIncThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to pbIncThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hController.stepThresholdVal(1);


function etFilterWindow_Callback(hObject, eventdata, handles)
% hObject    handle to etFilterWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etFilterWindow as text
%        str2double(get(hObject,'String')) returns contents of etFilterWindow as a double
try
  handles.hController.filterWindowSpec = str2num(get(hObject,'String'));
catch ME
  handles.hController.changeFilterSpec(); %Refresh view
end


% --- Executes during object creation, after setting all properties.
function etFilterWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etFilterWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in tblRareAdjustProps.
function tblRareAdjustProps_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tblRareAdjustProps (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes on button press in tbTab1.
function tbTab1_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab1
handles.hController.changeTabDisplayed(1);


% --- Executes on button press in tbTab2.
function tbTab2_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab2
handles.hController.changeTabDisplayed(2);

% --- Executes on button press in tbTab3.
function tbTab3_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab3
handles.hController.changeTabDisplayed(3);


% --- Executes on button press in tbTab4.
function tbTab4_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab4
handles.hController.changeTabDisplayed(4);


% --- Executes on button press in tbTab4.
function tbTab5_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab4
handles.hController.changeTabDisplayed(5);


% --- Executes on button press in tbTab4.
function tbTab6_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab4
handles.hController.changeTabDisplayed(6);


% --- Executes on button press in tbTab4.
function tbTab7_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab4
handles.hController.changeTabDisplayed(7);


% --- Executes on button press in tbTab4.
function tbTab8_Callback(hObject, eventdata, handles)
% hObject    handle to tbTab4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tbTab4
handles.hController.changeTabDisplayed(8);


% --- Executes on button press in tbShowAllProps.
function tbShowAllProps_Callback(hObject, eventdata, handles)
% hObject    handle to tbShowAllProps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%most.gui.toggleAdvancedPanel(hObject,16,'y');
hFig = ancestor(hObject,'figure');
most.gui.AdvancedPanelToggler.toggle(hFig);

% --- Executes on selection change in pmStimEventTypeDisplayed.
function pmStimEventTypeDisplayed_Callback(hObject, eventdata, handles)
% hObject    handle to pmStimEventTypeDisplayed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmStimEventTypeDisplayed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmStimEventTypeDisplayed
handles.hController.changeStimEventTypeDisplayed(hObject);

% --- Executes during object creation, after setting all properties.
function pmStimEventTypeDisplayed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmStimEventTypeDisplayed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbWaveformDisplay.
function rbWaveformDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to rbWaveformDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbWaveformDisplay
handles.hModel.displayMode = 'waveform';


% --- Executes on button press in rbRasterDisplay.
function rbRasterDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to rbRasterDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbRasterDisplay
handles.hModel.displayMode = 'raster';


% --- Executes on selection change in pmSpikeDisplayUnits.
function pmSpikeDisplayUnits_Callback(hObject, eventdata, handles)
% hObject    handle to pmSpikeDisplayUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmSpikeDisplayUnits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmSpikeDisplayUnits
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function pmSpikeDisplayUnits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmSpikeDisplayUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pbStimClassify.
function pbStimClassify_Callback(hObject, eventdata, handles)
% hObject    handle to pbStimClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hController.pbStimClassifyPressed();

% --- Executes on selection change in pmThresholdUnits.
function pmThresholdUnits_Callback(hObject, eventdata, handles)
% hObject    handle to pmThresholdUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmThresholdUnits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmThresholdUnits
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function pmThresholdUnits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmThresholdUnits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function mnuLoadCfg_Callback(hObject, eventdata, handles)
% hObject    handle to mnuLoadCfg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hModel.loadConfig();

% --------------------------------------------------------------------
function mnuSaveCfgAs_Callback(hObject, eventdata, handles)
% hObject    handle to mnuSaveCfgAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hModel.saveConfigAs();


% --------------------------------------------------------------------
function mnuFileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to mnuFileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when entered data in editable cell(s) in tblFreqAdjustProps.
function tblFreqAdjustProps_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tblFreqAdjustProps (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes on button press in pbPlotPSTH.
function pbPlotPSTH_Callback(hObject, eventdata, handles)
% hObject    handle to pbPlotPSTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.hModel.plotPSTH();


function etTimeBinMs_Callback(hObject, eventdata, handles)
% hObject    handle to etTimeBinMs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etTimeBinMs as text
%        str2double(get(hObject,'String')) returns contents of etTimeBinMs as a double
handles.hController.updateModel(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function etTimeBinMs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etTimeBinMs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function etPSTHAmpRange_Callback(hObject, eventdata, handles)
% hObject    handle to etPSTHAmpRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etPSTHAmpRange as text
%        str2double(get(hObject,'String')) returns contents of etPSTHAmpRange as a double
handles.hController.updateModel(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function etPSTHAmpRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etPSTHAmpRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
handles.hModel.quit();
delete(hObject);


% --------------------------------------------------------------------
function mnuViewMenu_Callback(hObject, eventdata, handles)
% hObject    handle to mnuViewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuViewWaveformDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to mnuViewWaveformDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnuViewRasterDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to mnuViewRasterDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
