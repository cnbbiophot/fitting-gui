function varargout = gui_FCS_SPA(varargin)
% GUI_FCS_SPA MATLAB code for gui_FCS_SPA.fig
%      GUI_FCS_SPA, by itself, creates a new GUI_FCS_SPA or raises the existing
%      singleton*.
%
%      H = GUI_FCS_SPA returns the handle to a new GUI_FCS_SPA or the handle to
%      the existing singleton*.
%
%      GUI_FCS_SPA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FCS_SPA.M with the given input arguments.
%
%      GUI_FCS_SPA('Property','Value',...) creates a new GUI_FCS_SPA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_FCS_SPA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_FCS_SPA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_FCS_SPA

% Last Modified by GUIDE v2.5 08-Apr-2019 14:33:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_FCS_SPA_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_FCS_SPA_OutputFcn, ...
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

% --- Executes just before gui_FCS_SPA is made visible.
function gui_FCS_SPA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_FCS_SPA (see VARARGIN)

% GUI_FCSfit = findobj(allchild(0), 'gui_FCSfit');
% handles_FCSfit = guidata(GUI_FCSfit);

if numel(varargin) == 6
    v.paramName = varargin{1};
    v.numParam = varargin{2};
    v.globalFit = varargin{3};
    v.dataSetSelection = varargin{4}; % Used as the matrix of all the datasets
    v.parFixedT = varargin{5};
    v.parSharedT = varargin{6};
    % Initialize variables
    v.dataset_selection = nan; % Initialize v.dataset_selection and v.searchRange 
    v.significance = str2num(get (handles.edit_significance, 'String'));
    v.searchRange = 10; 
    v.parameter_selection = nan; % When there has not been a selection yet
    % Set boxes
    set (handles.edit_searchRange, 'String', num2str(v.searchRange));
    set (handles.edit_CI, 'String', '')
    set (handles.pushbutton_plotCI, 'Enable', 'off')
elseif numel(varargin) > 6
    v.paramName = varargin{1};
    v.numParam = varargin{2};
    v.globalFit = varargin{3};
    v.dataSetSelection = varargin{4}; % Used as the pataset selected
    v.parFixedT = varargin{5};
    v.parSharedT = varargin{6};
    v.paramUnits = varargin{7};
    v.confidenceInterval = varargin{8};
    v.data = varargin{9}; % Data to plot the F-statistic
    v.Fchi_critical = varargin{10};
    v.parameter_selection = varargin{11};
    v.selectedDataset = varargin{12};
    v.sel_significance = varargin{13};
    v.sel_searchRange = varargin{14};
    % Initialize variables with former values
    v.dataset_selection = v.selectedDataset;
    v.significance = v.sel_significance;
    v.searchRange = v.sel_searchRange;
    % Fill the boxes with data
    % Disable boxes
    set (handles.edit_CI, 'String', num2str(v.confidenceInterval))
    set (handles.edit_significance, 'String', num2str(v.significance))
    set (handles.edit_searchRange, 'String', num2str(v.searchRange))
    set (handles.edit_significance, 'Enable', 'off')
    set (handles.pushbutton_SPA, 'Enable', 'off')
    set (handles.edit_searchRange, 'Enable', 'off')
end

setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

preparaParamTableVisible(v.globalFit, handles.table_selectParameters);
v.paramTable_allDataSets = inicializaParamTable(v.paramName, v.numParam, v.globalFit, v.dataSetSelection);
actualizaParamTableVisible(v.paramTable_allDataSets, v.parFixedT, v.parSharedT, v.paramName, v.numParam, v.globalFit, v.dataSetSelection, handles.table_selectParameters, v.parameter_selection, handles)

% Choose default command line output for gui_FCS_SPA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_FCS_SPA wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_FCS_SPA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
v=getappdata (handles.figure1, 'v'); %Recupera variables

varargout{1} = v.dataset_selection;
varargout{2} = v.parameter_selection;
varargout{3} = v.significance;
varargout{4} = v.searchRange;

delete(handles.figure1) % It closes the window



function edit_significance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_significance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_significance as text
%        str2double(get(hObject,'String')) returns contents of edit_significance as a double

v=getappdata (handles.figure1, 'v'); %Recupera variables

sign_var = str2double(get (handles.edit_significance, 'String'));

if or (sign_var > 1, sign_var < 0)
    disp('Significance cannot be greater than 1 or smaller than 0')
    v.significance = 0.32;
else
    v.significance = sign_var;
end    

setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables



% --- Executes during object creation, after setting all properties.
function edit_significance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_significance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CI as text
%        str2double(get(hObject,'String')) returns contents of edit_CI as a double


% --- Executes during object creation, after setting all properties.
function edit_CI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_SPA.
function pushbutton_SPA_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to pushbutton_SPA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v=getappdata (handles.figure1, 'v'); %Recupera variables

paramTable = get (handles.table_selectParameters, 'Data');   % Recover data from table

switch v.globalFit
    case false
        paramSelected=cell2mat(paramTable(:, 3));
        positionSel = find(paramSelected==1);
        if numel(positionSel) > 1 % Check if there are more than one chosen
            disp('Only one parameter can be chosen to fit, we choose the first on the list')
            positionSel = positionSel(1);
        end
        v.dataset_selection = 1;
        v.parameter_selection = positionSel;     
    case true
        paramSelected = cell2mat(paramTable(:, 5));
        positionSel = find(paramSelected==1);
        if numel(positionSel) > 1 % Check if there are more than one chosen
            disp('Only one parameter can be chosen to fit, we choose the first on the list')
            positionSel = positionSel(1);
        end
        
        if isempty(positionSel)
            v.dataset_selection = nan;
            v.parameter_selection = nan;
        else
            v.dataset_selection = paramTable{positionSel, 1}; % Set the output
            v.parameter_selection = positionSel;
        end
end


setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables
uiresume(handles.figure1) % Continues executing with new vargout

function infoColumns=preparaParamTableVisible(globalFit, h_tabla)
%Prepara la tabla con títulos y la limpia antes de empezar a rellenarla
%infoColumns=get (h_organelleTableGUI.uitable_organelleAnalysis, 'ColumnName'); %La primera fila contiene los nombres de las columnas. En el futuro la segunda podrá contener las unidades, comentarios, etc

if globalFit
    infoColumns={'Dataset', 'Parameter',...
         'Shared', 'Fixed', 'Select parameter'};
    columnFormat={'numeric', 'char'...
        'logical', 'logical', 'logical'};
    columnEditable= logical([0, 0, 0, 0, 1]);
else
    infoColumns={'Parameter', 'Fixed', 'Select parameter'};
    columnFormat={'char', 'logical', 'logical'};
    columnEditable= logical([0, 0, 1]);
end
set  (h_tabla, 'ColumnName', infoColumns)
set  (h_tabla, 'ColumnFormat', columnFormat)
set  (h_tabla, 'ColumnEditable', columnEditable)
set  (h_tabla, 'Data', [])

function paramTable_allDataSets=inicializaParamTable(paramName, numParam, globalFit, dataSetSelection)

% Inicializa la variable con los parámetros de la tabla
    idxdatasets=find(dataSetSelection);
    numDataSetsAjuste=numel(idxdatasets);

switch globalFit
    case false
    paramTable_allDataSets=cell(1, 3); %Si no hay globalFit
    filasParam=1:numParam;
    paramFijo=false(numParam ,1);
    paramSelected=false(numParam ,1);
    for n=1:numParam
        paramTable_allDataSets(filasParam(n), :)={paramName{n}, paramFijo(n), paramSelected(n)};
    end
    case true
        paramTable_allDataSets=cell(numDataSetsAjuste, 5); %Si hay globalFit
    for dataAjuste=1:numDataSetsAjuste
        filasParam=1+(dataAjuste-1)*numParam:(dataAjuste*numParam);
        paramFijo=false(numParam ,1);
        paramShared=false(numParam ,1);
        paramSelected=false(numParam ,1);
        for n=1:numParam
            paramTable_allDataSets(filasParam(n), :)={dataAjuste, paramName{n}, paramShared(n), paramFijo(n), paramSelected(n)};
        end
    end    
end

function actualizaParamTableVisible(paramTable_allDataSets, parFixedT, parSharedT, paramName, numParam, globalFit, dataSetSelection, h_tabla, parameter_selection, handles)
%Rellena la tabla con las filas de los parámetros de los datasets que
%mostrará


    idxdatasets=find(dataSetSelection);
    numDataSetsAjuste=numel(idxdatasets);

switch globalFit
    case false   
    paramTable=cell(1, 3); 
        filasParam=1:numParam; %Las filas que contienen los parámetros del ajuste
        paramFijo=cell2mat(parFixedT);
        paramSelected=cell2mat(paramTable_allDataSets(filasParam, 3));
        if ~isnan(parameter_selection) % check if parameter selection is available (second call), in that case, put it as such
            parameter_selection_dum = false(numParam,1);
            parameter_selection_dum(parameter_selection) = true;
            paramSelected = parameter_selection_dum(filasParam);
            clear parameter_selection_dum 
        end
        filasParam_new=1:numParam; %Las filas que contienen los parámetros del ajuste
        for n=1:numParam
            paramTable(filasParam_new(n), :)={paramName{n}, paramFijo(n), paramSelected(n)};
        end

    case true
        
    paramTable=cell(numDataSetsAjuste, 5); 
    paramShared_dum  = cell2mat(parSharedT);
    paramFijo_dum = cell2mat(parFixedT);
    for dataAjuste=1:numDataSetsAjuste
        filasParam = 1+(idxdatasets(dataAjuste)-1)*numParam:idxdatasets(dataAjuste)*numParam; %Las filas que contienen los parámetros del ajuste
        paramSelected = cell2mat(paramTable_allDataSets(filasParam, 5));
        if ~isnan(parameter_selection) % check if parameter selection is available (second call), in that case, put it as such
            parameter_selection_dum = false(numDataSetsAjuste*numParam,1);
            parameter_selection_dum(parameter_selection) = true;
            paramSelected = parameter_selection_dum(filasParam);
            clear parameter_selection_dum 
        end
        paramShared = paramShared_dum(filasParam);
        paramFijo = paramFijo_dum(filasParam);
        filasParam_new = 1+(dataAjuste-1)*numParam:(dataAjuste*(numParam)); %Las filas que contienen los parámetros del ajuste
        for n=1:numParam
            paramTable(filasParam_new(n), :)={paramTable_allDataSets{filasParam(n), 1}, paramName{n},...
                paramShared(n), paramFijo(n), paramSelected(n)};
        end
    end
end

set (h_tabla, 'Data', paramTable)

function  plot_CI(data, Fchi_critical, paramName, paramUnit, handles)

v=getappdata (handles.figure1, 'v'); %Recupera variables

h_axes = findobj('Tag', 'axes_SPA');
% set (h_axes, 'xlabel', [paramName ' / ' paramUnit]);
% set (h_axes, 'ylabel', '$\Chi^2/\Chi^2_{min}$', 'interpreter','latex');

h_axes = plot(data{1},data{2},'g', 'LineWidth', 3); % Plot data
hold on
h_axes = plot([min(data{1}) max(data{1})],[Fchi_critical Fchi_critical],'r', 'LineWidth', 3); % Plot line at Fchi_critical
hold off

pos_rel_selected_data = v.numParam - (v.numParam*v.selectedDataset - v.parameter_selection);
xlabel([paramName{pos_rel_selected_data} ' / (' paramUnit{pos_rel_selected_data} ')']);
ylabel( '$\chi^2/\chi^2_{min}$', 'interpreter','latex');

% --- Executes during object creation, after setting all properties.
function axes_SPA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_SPA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_SPA

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v=getappdata (handles.figure1, 'v'); %Recupera variables

if numel(varargin) > 6
    v.parameter_selection = nan;
    v.significance = nan;
elseif isnan(v.dataset_selection)
    v.parameter_selection = nan;
    v.significance = nan;
    disp('No parameter selected. Exiting Supported Plane Analysis')
end

setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

uiresume(handles.figure1)

function edit_searchRange_Callback(hObject, eventdata, handles)
% hObject    handle to edit_searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_searchRange as text
%        str2double(get(hObject,'String')) returns contents of edit_searchRange as a double

v=getappdata (handles.figure1, 'v'); %Recupera variables
searchRange = str2double(get (handles.edit_searchRange, 'String'));

v.searchRange = searchRange;
setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

% --- Executes during object creation, after setting all properties.
function edit_searchRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_searchRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_plotCI.
function pushbutton_plotCI_Callback(hObject, eventdata, handles, varargin)
% hObject    handle to pushbutton_plotCI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v=getappdata (handles.figure1, 'v'); %Recupera variables
plot_CI(v.data, v.Fchi_critical, v.paramName, v.paramUnits, handles);



% --- Executes during object creation, after setting all properties.
function pushbutton_plotCI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotCI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
