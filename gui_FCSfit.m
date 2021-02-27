function varargout = gui_FCSfit(varargin)
% [allParam chi2 dataSetSelection fittingFunction Gmodel]=gui_FCSfit(expData, fittingFunction, h_fitAxes, h_resAxes)
% allParam: todos los parámetros
% chi2: el chi2
% dataSetSelection son los índices de los datasets del ajuste (cuando hay más de uno)
% fittingFunction es la función de entrada
% Gmodel es la ACF para los datasets ajustados
%
% expData es una celda con los datos experimentales:
%   expData{1}={xdata1, ydata1, errdata1}
%   expData{N}={xdataN, ydataN, errdataN}
% fittingFunction es la función (por ahora sólo string)
% h_fitAxes, h_resAxes son los handles a los ejes de las gráficas del
% ajuste y los residuos (puede estar vacío)
%
%
% jri - 12-Feb-2015
% jri - 10-Oct-2018
% agv - 20-Nov-2018
% agv - 19Dic2018 Chi2 es correcto pero los errores asintóticos no, set them as Nan
% agv - 04Mar2019 Make function listbox functional
% agv - 14Mar2019 Implement copy buttons
% agv - 17jun2020 added log scale in Y axes

% % Example of use:
% % load('D:\Usuarios\Arturo\180307 FCS Volume Calibration\B&H data\analisis int10freq1000sin18\180307 633_10nm+t5_m1.mat')
% % expData={Gintervalos(:,:,:)};
% % gui_FCSfit(expData, 'fitfcn_FCS_3DTauD', [], [])
% % unicamente metiendo los datos
% % gui_FCSfit(expData)

% Last Modified by GUIDE v2.5 10-Feb-2015 15:08:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_FCSfit_OpeningFcn, ...
    'gui_OutputFcn',  @gui_FCSfit_OutputFcn, ...
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


% --- Executes just before gui_FCSfit is made visible.
function gui_FCSfit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_FCSfit (see VARARGIN)

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
ScreenUnits=get(0,'Units');
set(0,'Units','pixels');
ScreenSize=get(0,'ScreenSize');
set(0,'Units',ScreenUnits);
if isempty(gcbf)
    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
        GCBFPos(2) - FigHeight];
end

%set (handles.radiobutton_globalFit, 'Enable', 'off')
%set (handles.pushbutton_done, 'Enable', 'off')

v.globalFit=false;
v.numDataSets=numel(varargin{1});
v.data=varargin{1};

%v.data = get_uigetfile(); % gets the variable with a opening window

if numel(varargin) == 4
    v.fittingFunctionName=varargin{2};
    v.h_fitAxes=varargin{3}; %Handle a los ejes del ajuste
    v.h_resAxes=varargin{4}; %Handle a los ejes de los residuos
elseif numel(varargin) == 1
    v.fittingFunctionName='fitfcn_FCS_3D_G0_D_single';
    v.h_fitAxes=[]; %Handle a los ejes del ajuste
    v.h_resAxes=[]; %Handle a los ejes de los residuos
else
    error('Not enough input arguments')
end

v.dataSetSelection=true(v.numDataSets,1);
v.int_ajuste=zeros(v.numDataSets, 2);
for n=1:v.numDataSets
    v.int_ajuste(n, 1)=1;
    v.int_ajuste(n, 2)=size(v.data{n}, 1);
end

set (handles.edit_numDataSets, 'String', num2str(v.numDataSets))
set (handles.edit_chi2, 'String', '')

[v.paramName, v.paramUnits, v.numParam, v.fittingFunction]=parsefittingfunctionfile(v.fittingFunctionName);
set (handles.edit_fittingFunction, 'String', v.fittingFunctionName)
preparaParamTableVisible(v.globalFit, handles.table_fitParameters);
v.paramTable_allDataSets=inicializaParamTable(v.paramName, v.paramUnits, v.numParam, v.globalFit, v.numDataSets);
actualizaParamTableVisible(v.paramTable_allDataSets, v.paramName, v.paramUnits, v.numParam, v.globalFit, v.dataSetSelection, handles.table_fitParameters, handles);

v.allParam=cell(v.numDataSets, 1);
v.Gmodel=cell(v.numDataSets, 1);
v.chi2=zeros(v.numDataSets, 1);


v.colores=linspecer(v.numDataSets);
%Inicializa plots
v.h_fitPlot=zeros(v.numDataSets, 1);
v.h_resPlot=zeros(v.numDataSets, 1);
v.h_dataPlot=zeros(v.numDataSets, 1);
% La siguiente funcion inicializa los plots
[v.h_fig, v.h_fitAxes, v.h_resAxes, v.h_dataPlot, v.h_fitPlot, v.h_resPlot]=inicializaPlots (v.data, v.dataSetSelection, v.colores, v.h_fitAxes, v.h_resAxes, v.h_dataPlot, v.h_fitPlot, v.h_resPlot);

% v.colores(1,:)=[1 131 95]/255; %Verde
% v.colores(2,:)= [197 22 56]/255; %Rojo
% v.colores(3,:)= [0 102 204]/255; %Azul
% v.colores(4,:)= [50 50 50]/255; %Negro

setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

% Choose default command line output for gui_FCSfit
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_FCSfit wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_FCSfit_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
v=getappdata (handles.figure1, 'v'); %Recupera variables

varargout{1} = v.allParam;
varargout{2} = v.chi2;
varargout{3} = v.dataSetSelection;
varargout{4} = v.fittingFunction;
varargout{5} = v.Gmodel;
varargout{6} = handles.output;
setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables
delete(hObject);

function edit_numDataSets_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numDataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numDataSets as text
%        str2double(get(hObject,'String')) returns contents of edit_numDataSets as a double


% --- Executes during object creation, after setting all properties.
function edit_numDataSets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numDataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_fittingFunction_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fittingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function listbox_fitFun_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fittingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read box and selected item
listBoxItems = handles.listbox_fitFun.String;
listBoxSelectedIndexes = handles.listbox_fitFun.Value;
selectedString = listBoxItems{listBoxSelectedIndexes};

% get appdata
v=getappdata (handles.figure1, 'v');
v.fittingFunctionName = selectedString;

% write selected item in fittingfunction name
set (handles.edit_fittingFunction, 'String', selectedString)
% push button fittingfunction automatically
pushbutton_fittingFunction_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function edit_fittingFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fittingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_fit.
function pushbutton_fit_Callback(hObject, eventdata, handles)
% handles    structure with handles and user data (see GUIDATA)
v=getappdata (handles.figure1, 'v'); % Recupera variables (definidas en OpneningFunc)

tic

set (handles.pushbutton_done, 'Enable', 'on')   % Set the GUI to wait
set (handles.figure1,'Pointer','watch')
drawnow update

if handles.radiobutton_FCS.Value
    set (v.h_fitAxes, 'Yscale', 'linear') % Set the Y scale logatithmic 
    set (v.h_resAxes, 'Yscale', 'linear') % Set the Y scale linear
    set (v.h_fitAxes, 'Xscale', 'log') % Set the Y scale logatithmic 
    set (v.h_resAxes, 'Xscale', 'log') % Set the Y scale linear
else
    set (v.h_fitAxes, 'Yscale', 'log') % Set the Y scale logatithmic 
    set (v.h_resAxes, 'Yscale', 'log') % Set the Y scale linear
    set (v.h_fitAxes, 'Xscale', 'linear') % Set the Y scale logatithmic 
    set (v.h_resAxes, 'Xscale', 'linear') % Set the Y scale linear
end

%FUN está definida en parsefittingfunctionfile(v.fittingFunctionName);

idxdatasets=find(v.dataSetSelection);   % Indices de cada data set (si hay mas de 1)
numDataSetsAjuste=numel(idxdatasets);   % Numero de indices
paramTable=get (handles.table_fitParameters, 'Data');   % Los datos de la tabla

switch v.globalFit
    case false  % CASE NO GLOBAL
    for dataAjuste=1:numDataSetsAjuste  % Se hace un bucle para pasar por cada set de datos (en el global no habria que hacerlo)

        Gdata=cell2mat(v.data(idxdatasets(dataAjuste)));    % Toma los datos del set idxdatasets(dataAjuste)
        xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
        ydata=Gdata(:,2);
        yerr=Gdata(:,3);
        
            if ~all(yerr) && any(isnan(yerr)) % If any is zero, then make them 1 so no problems arise
                disp('Error are zero or Nan values. SET THEM TO 1, BUT THE FIT MAY BE WRONG')
                yerr = ones(size(yerr));
            end
        
        filasAllParam=(1+(idxdatasets(dataAjuste)-1)*v.numParam):(idxdatasets(dataAjuste)*v.numParam); %Las filas que contienen los parámetros del ajuste en allParam (y en paramTableAllDataSets)
        filasParamTable=(1+(dataAjuste-1)*v.numParam):(dataAjuste*v.numParam); 
        %Estos parámetros los coge directamente de la tabla, por si han cambiado
        paramFijo=cell2mat(paramTable(filasParamTable, 5)); % Parametros que son fijos (en la col 5 están los ticks)
        paramLibre=not(paramFijo); %Esto son logical, recuerda. = ¬paramFijo
        valorparametro=cell2mat(paramTable(filasParamTable, 6));    % Toma los valores de los parametros en la tabla
        valorLB=cell2mat(paramTable(filasParamTable, 8));   % Lower bound value
        valorUB=cell2mat(paramTable(filasParamTable, 9));   % Upper bound value

        [fitParam, allParam, chi2, deltaAllParam, ymodel]=...
            fitcore (v.fittingFunction, v.numParam, xdata, ydata, yerr, paramLibre, valorparametro, valorLB, valorUB);

        % Por el momento, los errores están mal calculados: los quitamos
        deltaAllParam = nan(size(deltaAllParam));
        
        %Mete los datos en la tabla que contiene todos los valores (paramTable_allDataSets)
        for n=1:v.numParam
            v.paramTable_allDataSets(filasAllParam(n), :)={paramTable{filasParamTable(n), 1}, v.paramName{n}, '', v.paramUnits{n},...
                paramFijo(n), allParam(n), deltaAllParam(n), valorLB(n), valorUB(n)};
        end

        v.allParam(idxdatasets(dataAjuste))={allParam}; % Guarda los nuevos parámetros
        v.chi2(idxdatasets(dataAjuste))=chi2;           % Guarda el nuevo chi2
        v.Gmodel(idxdatasets(dataAjuste))={[xdata ymodel yerr]}; % Guarda el modelo con los nuevos parámetros

        % Esto tiene que salir de aquí más tarde. Cuidado con ydata!!
        set (0, 'CurrentFigure', v.h_fig)

        if v.h_fitPlot(idxdatasets(dataAjuste))
            delete(v.h_fitPlot(idxdatasets(dataAjuste)));
            delete(v.h_resPlot(idxdatasets(dataAjuste)));
            v.h_fitPlot(idxdatasets(dataAjuste))=0;
            v.h_resPlot(idxdatasets(dataAjuste))=0;
        end
        v.h_fitPlot(idxdatasets(dataAjuste))=plot (v.h_fitAxes, xdata, ymodel, 'Color', v.colores(idxdatasets(dataAjuste), :), 'Linewidth', 2);
        v.h_resPlot(idxdatasets(dataAjuste))=plot (v.h_resAxes, xdata, ymodel-ydata, 'Color', v.colores(idxdatasets(dataAjuste), :), 'Linewidth', 2);
        ylim(v.h_resAxes, 'auto')

        set(handles.edit_chi2,'string',num2str(chi2)); % write the value of chi2 in the box

    end
    
    case true   % CASE GLOBAL 
        
        Gdata=v.data(idxdatasets)';    % Toma los datos del set idxdatasets(dataAjuste) ES UNA CELDA       
        % Toma los datos de las columnas
        paramShared=cell2mat(paramTable(:, 5));
        paramFijo=cell2mat(paramTable(:, 6)); % Parametros que son fijos (en la col 5 están los ticks)
        paramLibre=not(paramFijo);
        valorparametro=cell2mat(paramTable(:, 7));    % Toma los valores de los parametros en la tabla
        valorLB=cell2mat(paramTable(:, 9));   % Lower bound value
        valorUB=cell2mat(paramTable(:, 10));   % Upper bound value
        
        % Es el primer set de datos el que marca qué es global, se pone tick en el resto
        % Si una variable es global, su valor de frontera y si es fija o no
        % se repite en todos los dataset
        for n = 1:v.numParam % Recorre cada uno de los parámetros
            if cell2mat(paramTable(n, 5)) == 1 
                for k = 2:numel(idxdatasets) % Recorre los dataSet seleccionados a partir del 2º
                    fila_param = n +((k-1)*v.numParam);
                    paramTable(fila_param, 5) = {true};  % Tick global
                    paramTable(fila_param, 6) = paramTable(n, 6);    % Tick fijo
                    paramTable(fila_param, 9) = paramTable(n, 9);    % LB
                    paramTable(fila_param, 10) = paramTable(n, 10);  % UB
                end
            else
                for k = 2:numel(idxdatasets)
                    fila_param = n +((k-1)*v.numParam);
                    paramTable(n +((k-1)*v.numParam), 5) = {false}; 
                end
            end
        end
        
        % Actualizamos los parámetros tras hacer coincidir las variables
        % que fueran globales
        paramShared=cell2mat(paramTable(:, 5));
        paramFijo=cell2mat(paramTable(:, 6));
        paramLibre=not(paramFijo);
        valorLB=cell2mat(paramTable(:, 9));   
        valorUB=cell2mat(paramTable(:, 10));   

        [fitParam, ordered_allparam, allParam, chi2, deltaAllParam, ymodel]=...
            fitcore_global (v.fittingFunctionName, v.numParam, Gdata, paramShared, paramFijo, paramLibre, valorparametro, valorLB, valorUB);

        % Por el momento, los errores están mal calculados: los quitamos
        deltaAllParam = nan(size(deltaAllParam));
        
        % Mete los datos en la tabla que contiene todos los valores (paramTable_allDataSets)
        % Pero hay que reorganizar los datos, ya que salen del ajuste
        % global con otro orden
        paramName = repmat(v.paramName',numDataSetsAjuste,1); % In order to save the fitted variables
        paramUnits = repmat(v.paramUnits',numDataSetsAjuste,1);
        for n = 1:numel(idxdatasets) % Recorre cada uno de los set de datos
            for k = 1:v.numParam % Recorre cada variable
                fila_param_general = (v.numParam*(idxdatasets(n)-1) + 1) + (k-1); % Posición en la lista completa de parámetros
                f_p = (v.numParam*(n-1) + 1) + (k-1);   % Posición en la lista reducida de parámetros
                v.paramTable_allDataSets(fila_param_general, :)={paramTable{f_p, 1}, paramName{f_p}, '', paramUnits{f_p},...
                    paramShared(f_p), paramFijo(f_p), ordered_allparam(f_p), deltaAllParam(f_p), valorLB(f_p), valorUB(f_p)};
            end
        end
               
        for aj01=1:numel(Gdata) % Calcula el número de puntos en cada curva
            dataCurva=Gdata{aj01};
            numPtosxCurva(aj01)=size(dataCurva,1);
        end
        
        for dataAjuste=1:numDataSetsAjuste % Guarda los parámetros en la variable global y pinta las curvas ajustadas

            v.allParam(idxdatasets(dataAjuste))={ordered_allparam((1+(dataAjuste-1)*v.numParam:dataAjuste*v.numParam))}; % Guarda los nuevos parámetros 
             
            Gdata=cell2mat(v.data(idxdatasets(dataAjuste)));    % Toma los datos del set idxdatasets(dataAjuste)
            xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
            ydata=Gdata(:,2);
            yerr=Gdata(:,3);
            
                if ~all(yerr) && any(isnan(yerr)) % If any is zero, then make them 1 so no problems arise
                    disp('Error are zero values. SET THEM TO 1, BUT THE FIT MAY BE WRONG')
                    yerr = ones(size(yerr));
                end
                
            ymodel_set = ymodel(1+sum(numPtosxCurva(1:dataAjuste-1)):sum(numPtosxCurva(1:dataAjuste))); % Solo los puntos del dataset correspontiente
            v.Gmodel(idxdatasets(dataAjuste))={[xdata ymodel_set yerr]}; % Guarda el modelo con los nuevos parámetros

            % Esto tiene que salir de aquí más tarde. Cuidado con ydata!!
            set (0, 'CurrentFigure', v.h_fig)
                if v.h_fitPlot(idxdatasets(dataAjuste))
                    delete(v.h_fitPlot(idxdatasets(dataAjuste)));
                    delete(v.h_resPlot(idxdatasets(dataAjuste)));
                    v.h_fitPlot(idxdatasets(dataAjuste))=0;
                    v.h_resPlot(idxdatasets(dataAjuste))=0;
                end
                v.h_fitPlot(idxdatasets(dataAjuste))=plot (v.h_fitAxes, xdata, ymodel_set, 'Color', v.colores(idxdatasets(dataAjuste), :), 'Linewidth', 2);
                v.h_resPlot(idxdatasets(dataAjuste))=plot (v.h_resAxes, xdata, ymodel_set-ydata, 'Color', v.colores(idxdatasets(dataAjuste), :), 'Linewidth', 2);
                ylim(v.h_resAxes, 'auto')
        end
    v.chi2=chi2;           % Guarda el nuevo chi2
    set(handles.edit_chi2,'string',num2str(chi2));
end

%Actualiza la tabla de parámetros visible
actualizaParamTableVisible(v.paramTable_allDataSets, v.paramName, v.paramUnits, v.numParam, v.globalFit, v.dataSetSelection, handles.table_fitParameters, handles)

%set (handles.edit_chi2, 'String', num2str(v.chi2));

set (handles.figure1,'Pointer','arrow')
setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables
toc


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in pushbutton_copyTable.
function pushbutton_copyTable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=getappdata (handles.figure1, 'v');
PathName = uigetdir('D:\Usuarios\Arturo','Select folder to save data');

% Copy fit variables with names in a text file in rows in csv form
switch v.globalFit
    case false % CASE NO GLOBAL
        fit_data = v.paramTable_allDataSets(:,6);
        fit_data_table = cell(0,numel(v.paramName));
        for i = 1:v.numDataSets
            ind = 1 + (i-1) * v.numParam;
            fit_data_table = [fit_data_table;fit_data(ind:ind+v.numParam-1)']; 
        end
        fit_data = cell2table(fit_data_table, 'VariableNames',v.paramName);
        writetable((fit_data), [PathName '\FittedData.txt'])
    case true % CASE GLOBAL
        fit_data = v.paramTable_allDataSets(:,7);
        fit_data_table = cell(0,numel(v.paramName));
        for i = 1:v.numDataSets
            ind = 1 + (i-1) * v.numParam;
            fit_data_table = [fit_data_table;fit_data(ind:ind+v.numParam-1)']; 
        end
        fit_data = cell2table(fit_data_table, 'VariableNames',v.paramName);
        writetable((fit_data), [PathName '\FittedData.txt'])
end


% --- Executes on button press in pushbutton_copyCurves.
function pushbutton_copyCurves_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=getappdata (handles.figure1, 'v');
PathName = uigetdir('D:\Usuarios\Arturo','Select folder to save data');
exp_data2save = v.data;
fit_data2save = v.Gmodel;
% Create folder with date and the current time
t = [datetime('now', 'Format','yyMMdd-HHmmss')];
DateString = datestr(t,'yyMMdd-HHmmss');
folder_name = [PathName '\Fit_gui_FCSfit_' DateString];
mkdir(folder_name)
% Save as Matlab variable
save([folder_name '\data_fitting'],'exp_data2save', 'fit_data2save')
% Save as text
% Create different files for different datasets
idxdatasets=find(v.dataSetSelection);   % Indices de cada data set (si hay mas de 1)
numDataSetsAjuste=numel(idxdatasets);   % Numero de indices
for k = 1:numDataSetsAjuste    
    exp_data2save_dum = array2table([exp_data2save{k}(:,1)...
        exp_data2save{k}(:,2) exp_data2save{k}(:,3)],...
        'VariableNames',{'X','Y','Exp_Error'});
    fit_data2save_dum = array2table([fit_data2save{k}(:,1)...
        fit_data2save{k}(:,2) fit_data2save{k}(:,3)],...
        'VariableNames',{'X','Y','Exp_Error'});
    writetable((exp_data2save_dum), [folder_name '\data_experiment_' num2str(k) '.txt'])
    writetable((fit_data2save_dum), [folder_name '\data_fit_' num2str(k) '.txt'])
end

% --- Executes on button press in pushbutton_SupportPlaneAnalysis.
function pushbutton_SPA_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = getappdata (handles.figure1, 'v');
idxdatasets=find(v.dataSetSelection);   % Indices de cada data set (si hay mas de 1)

if v.globalFit % Select dataset to analyse if it is single fit
    % Do nothing
    parFixedT = v.paramTable_allDataSets(:,6);
    parSharedT = v.paramTable_allDataSets(:,5);
else % Select data from single dataset
    if v.numDataSets > 1
        dataSet_SPA = inputdlg('Select dataset:', 'Dataset for SPA', 1);
        dataSet_SPA = str2double(dataSet_SPA{1});
    else
        dataSet_SPA = 1;
    end
    filas_param = 1+(dataSet_SPA-1)*v.numParam : v.numParam + (dataSet_SPA-1)*v.numParam;
    parFixedT = v.paramTable_allDataSets(filas_param,5);
    parSharedT = [];
end

[dataset_selection, parameter_selection, sel_significance, sel_searchRange]...
    = gui_FCS_SPA(v.paramName, v.numParam, v.globalFit, v.dataSetSelection, parFixedT, parSharedT); % First call to the function

if isempty(parameter_selection)% Check if there are selected variables
    disp('No parameter selected')
    return
end

if v.globalFit % Select dataset to analyse if it is single fit
    valorparam_SPA = cell2mat(v.paramTable_allDataSets(parameter_selection, 7));
else
    valorparam_SPA = cell2mat(v.paramTable_allDataSets(parameter_selection, 6));
end  

searchRange_param = valorparam_SPA*sel_searchRange/100; % Search in of the value
disp(['Search range is from ' num2str(valorparam_SPA - searchRange_param) ' to ' num2str(valorparam_SPA + searchRange_param)])

% Compute confidence interval
[confidenceInterval, paramSupportPlane, Fchi, Fchi_critical]=...
    compute_SPA (parameter_selection, searchRange_param, sel_significance, handles, v, dataset_selection);

% Data to plot the Support-plane analysis
xSPA = paramSupportPlane; 
ySPA = Fchi;
dataSPA = {xSPA, ySPA};

[~,~,~,~] = gui_FCS_SPA(v.paramName, v.numParam, v.globalFit, v.dataSetSelection, ...
        parFixedT, parSharedT, v.paramUnits, confidenceInterval, dataSPA, Fchi_critical, ...
        parameter_selection, dataset_selection, sel_significance, sel_searchRange); % Second call to the function to plot

function edit_chi2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chi2 as text
%        str2double(get(hObject,'String')) returns contents of edit_chi2 as a double


% --- Executes during object creation, after setting all properties.
function edit_chi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chi2 (see GCBO)
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

%checkDoneState=get (handles.pushbutton_done, 'Enable');
%if strcmpi (checkDoneState, 'on')
uiresume(hObject);
%end

% --- Executes on button press in pushbutton_fittingFunction.
function pushbutton_fittingFunction_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_fittingFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=getappdata (handles.figure1, 'v'); %Recupera variables

v.fittingFunctionName=get (handles.edit_fittingFunction, 'String');
[v.paramName, v.paramUnits, v.numParam, v.fittingFunction]=parsefittingfunctionfile(v.fittingFunctionName);
set (handles.edit_fittingFunction, 'String', v.fittingFunctionName)

preparaParamTableVisible(v.globalFit, handles.table_fitParameters);
v.paramTable_allDataSets=inicializaParamTable(v.paramName, v.paramUnits, v.numParam, v.globalFit, v.numDataSets);
actualizaParamTableVisible(v.paramTable_allDataSets, v.paramName, v.paramUnits, v.numParam, v.globalFit, v.dataSetSelection, handles.table_fitParameters, handles);
v.allParam=cell(v.numDataSets, 1);
v.Gmodel=cell(v.numDataSets, 1);

setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

function infoColumns=preparaParamTableVisible(globalFit, h_tabla)
%Prepara la tabla con títulos y la limpia antes de empezar a rellenarla
%infoColumns=get (h_organelleTableGUI.uitable_organelleAnalysis, 'ColumnName'); %La primera fila contiene los nombres de las columnas. En el futuro la segunda podrá contener las unidades, comentarios, etc

if globalFit
    infoColumns={'Dataset', 'Parameter', 'Meaning', 'Units',...
        'Share', 'Fix', 'Value', 'Error', ...
        'Lower bound', 'Upper bound'};
    columnFormat={'numeric', 'char', 'char', 'char',...
        'logical', 'logical', 'numeric', 'numeric', ...
        'numeric', 'numeric'};
    columnEditable= logical([0, 0, 0, 0, 1, 1, 1, 0, 1, 1]);
    v.fittingFunctionName='fitfcn_FCS_3D_G0_D_z0_omega0';
else
    infoColumns={'Dataset', 'Parameter', 'Meaning', 'Units',...
        'Fix', 'Value', 'Error', ...
        'Lower bound', 'Upper bound'};
    columnFormat={'numeric', 'char', 'char', 'char',...
        'logical', 'numeric', 'numeric', ...
        'numeric', 'numeric'};
    columnEditable= logical([0, 0, 0, 0, 1, 1, 0, 1, 1]);
end
set  (h_tabla, 'ColumnName', infoColumns)
set  (h_tabla, 'ColumnFormat', columnFormat)
set  (h_tabla, 'ColumnEditable', columnEditable)
set (h_tabla, 'Data', [])

% --- Executes when selected object is changed in uipanel_fit.
function uipanel_fit_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_fit
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

v=getappdata (handles.figure1, 'v'); %Recupera variables
if eventdata.NewValue==handles.radiobutton_globalFit
    v.globalFit=true;
else
    v.globalFit=false;
end
preparaParamTableVisible(v.globalFit, handles.table_fitParameters);
setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables
% push button fittingfunction automatically
pushbutton_fittingFunction_Callback(hObject, eventdata, handles)

% --- Executes on button press in pushbutton_dataSets.
function pushbutton_dataSets_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=getappdata (handles.figure1, 'v'); %Recupera variables
[v.dataSetSelection, v.int_ajuste, v.numDataSets]=gui_FCSfit_choosedataset(v.data, v.int_ajuste, v.dataSetSelection, v.colores);
set (handles.edit_numDataSets, 'String', num2str(v.numDataSets));
actualizaParamTableVisible(v.paramTable_allDataSets, v.paramName, v.paramUnits, v.numParam, v.globalFit, v.dataSetSelection, handles.table_fitParameters,handles);
[v.h_fig, v.h_fitAxes, v.h_resAxes, v.h_dataPlot, v.h_fitPlot, v.h_resPlot]=inicializaPlots (v.data, v.dataSetSelection, v.colores, v.h_fitAxes, v.h_resAxes, v.h_dataPlot, v.h_fitPlot, v.h_resPlot);
setappdata(handles.figure1, 'v', v); %Guarda los cambios en variables

function [fitparam, allParam, chi2, deltaAllParam, ymodel]=fitcore (FUN, numParam, xdata, ydata, yerr, paramlibre, valorparametro, valorLB, valorUB)

indparamvariables=[];
indparamfijos=[];
valorparamfijos=[];
numparamfijos=0;
numparamvariables=0;
guess=[];
LB=[];
UB=[];

%Falta definir el intervalo de ajuste v.int_ajuste

for n=1:numParam
    [guess LB UB numparamvariables numparamfijos indparamvariables indparamfijos valorparamfijos]=...
        check_fijoovariable (paramlibre(n), valorparametro(n), n, guess, LB, UB, valorLB(n), valorUB(n), numparamvariables, numparamfijos, indparamvariables, indparamfijos, valorparamfijos);
end


[fitparam, allParam, resnorm, resfun, EXITFLAG, OUTPUT, LAMBDA, jacob_mat]...
    = ajusta_lsqnonlin(FUN, guess, LB ,UB ,[], xdata, ydata, yerr, indparamvariables, indparamfijos, valorparamfijos);

[chi2 deltaParamFit ymodel] = ajusta_computeuncertainties (FUN, xdata, ydata, yerr, allParam, indparamvariables, jacob_mat); %Simplemente calcula el modelo y las incertidumbres del ajuste

deltaAllParam=zeros(1, numParam);
deltaAllParam(paramlibre)=deltaParamFit;

function [fitparam, ordered_allparam, allParam, chi2, deltaAllParam, ymodel]=...
    fitcore_global (FUN, numParam, dataCell, varGlobals_table, varFix_table, paramlibre, paramsIni_table, valorLB_table, valorUB_table)
% La estructura de la función cambia ligeramente respecto al anterior
% fitcore, ya que está basada en ajusteGlobalGral
% La salida de esta función es respecto a los set de datos que se estén
% analizando

% Redimensionamos y recolocamos los vectores de variables para la función ajusteGral
num_datasets = numel(dataCell);
num_variables = num_datasets*numParam - (num_datasets-1)*numel(find(varGlobals_table==1));

varFix = zeros(num_variables,1);
paramsIni = zeros(num_variables,1);
valorLB = zeros(num_variables,1);
valorUB = zeros(num_variables,1); 

num_param_curr = 1;

for n = 1:numParam  % Se ponen las variables para que la utilice ajusteGlobalGral
    if varGlobals_table(n) == 1     % Si es global 
        varFix(num_param_curr) = varFix_table(n);
        paramsIni(num_param_curr) = paramsIni_table(n);
        valorLB(num_param_curr) = valorLB_table(n);
        valorUB(num_param_curr) = valorUB_table(n);
        num_param_curr = num_param_curr + 1;
    else    % Si no es global
        num_set = 1;
        for m = 1:num_datasets
            varFix(num_param_curr) = varFix_table(n+numParam*(m-1));
            paramsIni(num_param_curr) = paramsIni_table(n+numParam*(m-1));
            valorLB(num_param_curr) = valorLB_table(n+numParam*(m-1));
            valorUB(num_param_curr) = valorUB_table(n+numParam*(m-1));
            num_param_curr = num_param_curr + 1;
        end
    end
end

varGlobals = varGlobals_table(1:numParam);

[indparam, fitparam, allParam, indparamvariables, indparamfijos, jacob_mat, dsid]...
    = ajusteGlobalGral_gui(dataCell, varGlobals', paramsIni', varFix', valorLB', valorUB', FUN);

Gdata=cell2mat(dataCell);
xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
ydata=Gdata(:,2);
yerr=Gdata(:,3);

% % This was a test due to the fact that with batch correlate error were set
% % to zero
% %     if ~all(yerr) && any(isnan(yerr)) % If any is zero, then make them 1 so no problems arise
% %         disp('Error are zero values. SET THEM TO 1, BUT THE FIT MAY BE WRONG')
% %         yerr = ones(size(yerr));
% %     end

[chi2, deltaParamFit, ymodel] = ajusta_computeuncertainties_global (FUN, dsid, xdata, ydata, yerr, indparamvariables, indparamfijos, fitparam, allParam, indparam, jacob_mat); %Simplemente calcula el modelo y las incertidumbres del ajuste

% Se recolocan los parámetros ajustados para la salida
indices_old_params = reshape(indparam',1,numel(indparam));
deltaAllParam = zeros(1,numel(indparam));
for i = 1:numel(indparam)
    ordered_allparam(i) = allParam(indices_old_params(i));
    deltaAllParam(i) = abs(deltaParamFit(indices_old_params(i)));
end

function paramTable_allDataSets=inicializaParamTable(paramName, paramUnits, numParam, globalFit, numTotalDataSets)

% Inicializa la variable con los parámetros de la tabla

switch globalFit
    case false
    paramTable_allDataSets=cell(numTotalDataSets, 9); %Si no hay globalFit, creo que son sólo 10.
    for dataAjuste=1:numTotalDataSets
        filasParam=1+(dataAjuste-1)*numParam:(dataAjuste*numParam);
        paramFijo=false(numParam ,1);
        valorParametro=ones (numParam, 1);
        deltaParametro=zeros (numParam, 1);
        valorLB = zeros(numParam, 1);
        valorUB = Inf*ones(numParam, 1);
        for n=1:numParam
            paramTable_allDataSets(filasParam(n), :)={dataAjuste, paramName{n}, '', paramUnits{n}, paramFijo(n), valorParametro(n), deltaParametro(n), valorLB(n), valorUB(n)};
        end
    end
    case true
        paramTable_allDataSets=cell(numTotalDataSets, 10); %Si no hay globalFit, creo que son sólo 10.
    for dataAjuste=1:numTotalDataSets
        filasParam=1+(dataAjuste-1)*numParam:(dataAjuste*numParam);
        paramShared=false(numParam ,1);
        paramFijo=false(numParam ,1);
        valorParametro=ones (numParam, 1);
        deltaParametro=zeros (numParam, 1);
        valorLB = zeros(numParam, 1);
        valorUB = Inf*ones(numParam, 1);
        for n=1:numParam
            paramTable_allDataSets(filasParam(n), :)={dataAjuste, paramName{n}, '', paramUnits{n}, paramShared(n), paramFijo(n), valorParametro(n), deltaParametro(n), valorLB(n), valorUB(n)};
        end
    end    
end

function actualizaParamTableVisible(paramTable_allDataSets, paramName, paramUnits, numParam, globalFit, dataSetSelection, h_tabla, handles)
%Rellena la tabla con las filas de los parámetros de los datasets que
%mostrará

    idxdatasets=find(dataSetSelection);
    numDataSetsAjuste=numel(idxdatasets);

switch globalFit
    case false
    paramTable=cell(numDataSetsAjuste, 9); %Si no hay globalFit, creo que son sólo 9

    for dataAjuste=1:numDataSetsAjuste
        filasParam=1+(idxdatasets(dataAjuste)-1)*numParam:idxdatasets(dataAjuste)*numParam; %Las filas que contienen los parámetros del ajuste

        paramFijo=cell2mat(paramTable_allDataSets(filasParam, 5));
        valorParam=cell2mat(paramTable_allDataSets(filasParam, 6));
        deltaParam=cell2mat(paramTable_allDataSets(filasParam, 7));
        valorLB=cell2mat(paramTable_allDataSets(filasParam, 8));
        valorUB=cell2mat(paramTable_allDataSets(filasParam, 9));
        chi2=cell2mat(paramTable_allDataSets(filasParam, 7));
        filasParam_new=1+(dataAjuste-1)*numParam:(dataAjuste*(numParam)); %Las filas que contienen los parámetros del ajuste
        for n=1:numParam
            paramTable(filasParam_new(n), :)={paramTable_allDataSets{filasParam(n), 1}, paramName{n}, '', paramUnits{n},...
                paramFijo(n), valorParam(n), deltaParam(n), valorLB(n), valorUB(n)};
        end

    end
    
    case true
    paramTable=cell(numDataSetsAjuste, 10); % Con globalFit son 10
       
    for dataAjuste=1:numDataSetsAjuste
        filasParam=1+(idxdatasets(dataAjuste)-1)*numParam:idxdatasets(dataAjuste)*numParam; %Las filas que contienen los parámetros del ajuste

        paramShared  = cell2mat(paramTable_allDataSets(filasParam, 5));
        paramFijo=cell2mat(paramTable_allDataSets(filasParam, 6));
        valorParam=cell2mat(paramTable_allDataSets(filasParam, 7));
        deltaParam=cell2mat(paramTable_allDataSets(filasParam, 8));
        valorLB=cell2mat(paramTable_allDataSets(filasParam, 9));
        valorUB=cell2mat(paramTable_allDataSets(filasParam, 10));
        chi2=cell2mat(paramTable_allDataSets(filasParam, 6));
        filasParam_new=1+(dataAjuste-1)*numParam:(dataAjuste*(numParam)); %Las filas que contienen los parámetros del ajuste
        for n=1:numParam
            paramTable(filasParam_new(n), :)={paramTable_allDataSets{filasParam(n), 1}, paramName{n}, '', paramUnits{n},...
                paramShared(n), paramFijo(n), valorParam(n), deltaParam(n), valorLB(n), valorUB(n)};
        end
    end
end

set (h_tabla, 'Data', paramTable)

function [h_fig, h_fitAxes, h_resAxes, h_dataPlot, h_fitPlot, h_resPlot]=inicializaPlots (dataSet, dataSetSelection, colores, h_fitAxes_in, h_resAxes_in, h_dataPlot_in, h_fitPlot_in, h_resPlot_in)
% Borra o dibuja los puntos y las líneas de los datasets que están activados en cada momento
h_fitAxes=h_fitAxes_in;
h_resAxes=h_resAxes_in;
h_dataPlot=h_dataPlot_in;
h_fitPlot=h_fitPlot_in;
h_resPlot=h_resPlot_in;

idxdatasets=find(dataSetSelection);
numDataSetsAjuste=numel(idxdatasets);
numTotalDataSets=numel(dataSet);

if isempty(h_fitAxes)
    h_fig=figure;
    h_fitAxes=subplot (2, 1, 1);
    h_resAxes=subplot (2, 1, 2);
    set(h_fitAxes, 'position', [0.1,0.3,0.8,0.65] );
    set(h_resAxes, 'position', [0.1,0.1,0.8,0.1] );
else
    h_fig=get (h_fitAxes, 'Parent');
%    cla (h_fitAxes)
%    cla (h_resAxes)
end

hold (h_fitAxes, 'on')
hold (h_resAxes, 'on')

for n=1:numTotalDataSets
    if h_dataPlot(n)
        delete(h_dataPlot(n));
    end
    if and(h_fitPlot(n), not(dataSetSelection(n)))
        delete(h_fitPlot(n));
        delete(h_resPlot(n));
        h_fitPlot(n)=0;
        h_resPlot(n)=0;
    end
end

h_dataPlot=zeros(numTotalDataSets, 1); 
for dataAjuste=1:numDataSetsAjuste
    Gdata=cell2mat(dataSet(idxdatasets(dataAjuste)));
    xdata=Gdata(:,1);
    ydata=Gdata(:,2);
    yerr=Gdata(:,3);
    set(0, 'CurrentFigure', h_fig)
    set(h_fig, 'CurrentAxes', h_fitAxes)
    h_dataPlot(idxdatasets(dataAjuste))=errorbar(xdata, ydata, yerr, 'o', 'Color', colores(idxdatasets(dataAjuste), :), 'Linewidth', 2);
end

% set ([h_fitAxes, h_resAxes], 'Xscale', 'log')
set ([h_fitAxes], 'Xscale', 'linear') %added
set ([h_fitAxes], 'Yscale', 'log') %adde



function expData = get_uigetfile()
% Función para abrir un cuadro de diálogo para cargar la variable. Sólo
% carga una variable. Hay que modificarlo para poder añadir variables
    [FileName,PathName] = uigetfile();
    load([PathName FileName]);
    expData = Gintervalos(:,:,:);