function [confidenceInterval, paramSupportPlane, Fchi, Fchi_critical]=...
    compute_SPA (indParamSupport, searchRange, significance, handles, v, dataset_selection)

% Análsis del intervalo de confianza cambiando un parámetro de ajuste (indParamSupport)
% Este archivo está hecho pensando en un ajuste biexponencial, pero
% cambiando el modelo serviría para cualquier ajuste. 
% Se puede contrastar con el mismo que hice para Kd en 2016
%
% INPUT
% indParamSupport: is the index in the global variable vector of the
% variable to analyse in the Support-plane analysis
% searchRange: half width of the interval to search for the Support-plane
% analysis. Computed in percentage from gui_FCS_SPA
% significance: probability that the central value (fitting value) lies
% within the computed interval
% handles: handles inherited from gui_FCSfit (to retrieve table with data)
% v: appdata inherited from gui_FCSfit (to get the vectors with the
% variables fixed, global, etc)
% dataset_selection: number of dataset of the selected variable
%
% OUTPUT
% confidenceInterval: computed interval within which the real variable
% would lie with a percentage given by the significance (actual
% support-plane analysis)
% paramSupportPlane: whole interval to search in the analysed variable
% Fchi: plane of the F statistic (chi2SupportPlane/chi2_min) in the interval given by paramSupportPlane
% Fchi_critical: critical value in F statistic to search the interval SPA
%
% agv 12Abr2019

%Número de puntos en la abscisa (pármetro de soporte para el cálculo del intervalo de confianza)
% numPuntosParamSupport=500; % Number of points in the search interval
numPuntosParamSupport=50; % This yields reasonably accurate values

% Get optimal value of the parameter, number of points and the variable
% parameters
paramTable=get (handles.table_fitParameters, 'Data');   % Los datos de la tabla
switch v.globalFit
    case true
        valorparametro=cell2mat(paramTable(indParamSupport, 7));    % Toma los valores de los parametros en la tabla
        idxdatasets=find(v.dataSetSelection);   % Indices de cada data set (si hay mas de 1)
        Gdata=cell2mat(v.data(idxdatasets));    % Toma los datos del set idxdatasets(dataAjuste)
        xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
        paramFix=cell2mat(paramTable(:, 6));
        eta1=numel(find(paramFix==false)); % eta1 is the number of variables in the fit
        eta2=numel(xdata); % eta2 is the number of experimental points to fit
        chi2_min = v.chi2;
    case false
        valorparametro=cell2mat(paramTable(indParamSupport, 6));    % Toma los valores de los parametros en la tabla
        idxdatasets=find(dataset_selection);   % Indices de cada data set (si hay mas de 1)
        Gdata=cell2mat(v.data(idxdatasets));   % Toma los datos del set idxdatasets(dataAjuste)
        xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
        paramFix=cell2mat(paramTable(:, 5));
        eta1=numel(find(paramFix==false)); % eta1 is the number of variables in the fit
        eta2=numel(xdata); % eta2 is the number of experimental points to fit
        chi2_min = v.chi2(dataset_selection);
end

paramSupport_Fit = valorparametro; % Fitted value of the parameter to fit
limiteSuperior_paramSupport = paramSupport_Fit + searchRange;  % Limits of the supported parameter plane (paramSupportPlane)
limiteInferior_paramSupport = paramSupport_Fit - searchRange;
if limiteInferior_paramSupport<=0 % If the interval goes beyond zero, reescale it to end in zero
    limiteInferior_paramSupport=1E-15; % Para evitar problemas al resolver la ecuación
end
paramSupportPlane=linspace(limiteInferior_paramSupport, limiteSuperior_paramSupport, numPuntosParamSupport);

% Compute Fchi (Fchi=chi2SupportPlane/chi2_min) por every point in paramSupportPlane
[Fchi, ~]=Fchi_fcn(v, handles, dataset_selection, indParamSupport, paramSupportPlane, chi2_min);

% Look for minimum of the Fchi distribution and check whether it is in the
% center or not
[~, idx]=min(Fchi);
paramSupportMin=paramSupportPlane(idx);
if ~(paramSupportMin == paramSupport_Fit) % Check if it went correctly
    disp('Minimum of Support-plane differs from central value')
end

% Compute the Fchi critical thoretically (maximum value that takes the F
% distribution for a certain significance)
Fchi_critical = 1 + f_inv(1-significance, eta1,eta2)*eta1/eta2; % Sale de la pag 135 de Lakowic

% Resuelve para el parámetro (la incógnita es el parámetro , por eso param_x)
% Es decir, calculo el valor del parámetro para el cual función Fchi=chi2(par)/chi2min-Fchi_critical=0

% Set the equation (y=Fchi-Fchi_critical) to find the roots (the zero values) 
fun_paramSupportCI = @(param_x) paramSupportCIsolve(param_x, Fchi_critical, v, handles, ...
    dataset_selection, indParamSupport, chi2_min);    % function of param_x alone
%Comprueba que hay un cero (por cambio de signo). Si no, devuelve NaN
SupportPlane_Y = fun_paramSupportCI(paramSupportPlane);

% Check if there are zeros and compute the zeros
if or(sign(SupportPlane_Y(1))==sign(fun_paramSupportCI(paramSupport_Fit)) ,...
        sign(SupportPlane_Y(length(SupportPlane_Y)))==sign(fun_paramSupportCI(paramSupport_Fit)))
    paramSupportCI(1)=NaN;
    paramSupportCI(2)=NaN;
    s=sprintf('Chi2(par) no tiene mínimo en el intervalo [%f, %f] o el intervalo de confianza para alfa=%f es mayor',...
        limiteInferior_paramSupport, limiteSuperior_paramSupport, significance);
    disp (s)
else
    paramSupportCI(1) = fzero(fun_paramSupportCI, [limiteInferior_paramSupport, paramSupport_Fit]); %El intervalo de búsqueda es tres veces sigma, ya que mu±3*sigma contiene el 99.9% de la probabilidad de encontrar mu
    paramSupportCI(2) = fzero(fun_paramSupportCI, [paramSupport_Fit, limiteSuperior_paramSupport]);
    
end

confidenceInterval=paramSupportCI(2)-paramSupportCI(1);

%{
%Resuelvo numéricamente la intersección de la función Fchi con el valor de Fchi_critical
% Esto es encontrar los cambios de signo en Fchi
% Un problema que puede ocurrir es que no encuentre los cambios de signo debido a que el intervalo esté mal.
% Es muy lento porque para tener buena resolución tengo que evaluar Fchi en muchos puntos
temp=[0 diff(sign(Fchi-Fchi_critical))]; %Pongo un cero al principio porque al hacer diff se pierde un punto
zeroCrossings=find(temp);
KdConfidence(1)=KdSupportPlane(zeroCrossings(1));
KdConfidence(2)=KdSupportPlane(zeroCrossings(2));
confidenceInterval=KdConfidence(2)-KdConfidence(1);
%}

function [Fchi chi2SupportPlane]=Fchi_fcn(v, handles, dataset_selection, indParamSupport, paramSupportPlane, chi2_min)
% Función Fchi. Ed. 4.24 Lakowicz viejo:
% Fchi=chi2(par)/chi2(min)

% Allocate space
chi2SupportPlane=zeros(size(paramSupportPlane));

% Compute the plane of Chi2 within paramSupportPlane
for n=1:numel(paramSupportPlane)
    chi2SupportPlane(n) = fit_SPA(v, handles, dataset_selection, indParamSupport, paramSupportPlane(n));
end

% F distribution definition
Fchi=chi2SupportPlane/chi2_min;


function y=paramSupportCIsolve(paramSupportPlane, Fchi_critical, v, handles, dataset_selection, indParamSupport, chi2_min)
% Ecuación f(paramSupport)=Fchi(paramSupport)-Fchicritical=0 que uso para resolver para el parámetro para la que Fchi(paramSupport)=Fcritical
% Para eso tengo que calcular chi2(par) para cada valor del parámetro que varío (paramSupport), que se encuentra en el vector paramSupportPlane. Si fuese un ajuste
% de varios parámetros tendría que ajustar dejando tao fijo y variando el resto

%Calculo Fchi para el valor del parámetro correspondiente
[Fchi, ~]=Fchi_fcn(v, handles, dataset_selection, indParamSupport, paramSupportPlane, chi2_min);

y=Fchi-Fchi_critical; %Ésta es la función que minimiza


function chi2 = fit_SPA(v, handles, dataset_selection, indParamSupport, valueParamSupport) 
% Compute Chi2

idxdatasets=find(v.dataSetSelection);   % Indices de cada data set (si hay mas de 1)
paramTable=get (handles.table_fitParameters, 'Data');   % Los datos de la tabla

switch v.globalFit
        case false  % CASE NO GLOBAL
        dataAjuste = dataset_selection;  % Se hace un bucle para pasar por cada set de datos (en el global no habria que hacerlo)

        Gdata=cell2mat(v.data(idxdatasets(dataAjuste)));    % Toma los datos del set idxdatasets(dataAjuste)
        xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
        ydata=Gdata(:,2);
        yerr=Gdata(:,3);
        
            if ~all(yerr) && any(isnan(yerr)) % If any is zero, then make them 1 so no problems arise
                disp('Error are zero or Nan values. SET THEM TO 1, BUT THE FIT MAY BE WRONG')
                yerr = ones(size(yerr));
            end
        filasParamTable=(1+(dataAjuste-1)*v.numParam):(dataAjuste*v.numParam); 
        %Estos parámetros los coge directamente de la tabla, por si han cambiado
        paramFijo=cell2mat(paramTable(filasParamTable, 5)); % Parametros que son fijos (en la col 5 están los ticks)
        paramLibre=not(paramFijo); %Esto son logical, recuerda. = ¬paramFijo
        valorparametro=cell2mat(paramTable(filasParamTable, 6));    % Toma los valores de los parametros en la tabla
        
        % Make the evaluation parameter fix and fix the value
        paramLibre(indParamSupport) = false; % Make it fix (single fit works with paramLibre!)
        valorparametro(indParamSupport) = valueParamSupport; % Give the desired value
        
        valorLB=cell2mat(paramTable(filasParamTable, 8));   % Lower bound value
        valorUB=cell2mat(paramTable(filasParamTable, 9));   % Upper bound value

        chi2 = fitcore_SPA (v.fittingFunction, v.numParam, xdata, ydata, yerr, paramLibre, valorparametro, valorLB, valorUB);
        
        case true   % CASE GLOBAL 
        
        % Take all the variables from the table inherited from gui_FCSfit
        Gdata=v.data(idxdatasets)';    % Toma los datos del set idxdatasets(dataAjuste) ES UNA CELDA       
        % Toma los datos de las columnas
        paramShared=cell2mat(paramTable(:, 5));
        paramFijo=cell2mat(paramTable(:, 6)); % Parametros que son fijos (en la col 5 están los ticks)
        paramLibre=not(paramFijo);
        valorparametro=cell2mat(paramTable(:, 7));    % Toma los valores de los parametros en la tabla
        valorLB=cell2mat(paramTable(:, 9));   % Lower bound value
        valorUB=cell2mat(paramTable(:, 10));   % Upper bound value
        
        if paramShared(indParamSupport) == 1 % If the parameter to analyse is global, fix it in every dataset
            pos_rel_selected_data = v.numParam - (v.numParam*dataset_selection - indParamSupport);
            for d = 1:numel(idxdatasets)
                pos_abs = v.numParam*(d-1) + pos_rel_selected_data;
                paramFijo(pos_abs) = true;
                valorparametro(pos_abs) = valueParamSupport;
            end
        else % If it is not global
            paramFijo(indParamSupport) = 1; % Fix the parameter
            valorparametro(indParamSupport) = valueParamSupport; % Give the desired value
        end
        
        % Fit the equation again with the chosen parameter fixed and with
        % the desired value
        chi2 = fitcore_global_SPA (v.fittingFunctionName, v.numParam, Gdata, paramShared, paramFijo, paramLibre, valorparametro, valorLB, valorUB);

end
        
function chi2 = fitcore_SPA (FUN, numParam, xdata, ydata, yerr, paramlibre, valorparametro, valorLB, valorUB)

indparamvariables=[];
indparamfijos=[];
valorparamfijos=[];
numparamfijos=0;
numparamvariables=0;
guess=[];
LB=[];
UB=[];

for n=1:numParam
    [guess, LB, UB, numparamvariables, numparamfijos, indparamvariables, indparamfijos, valorparamfijos]=...
        check_fijoovariable (paramlibre(n), valorparametro(n), n, guess, LB, UB, valorLB(n), valorUB(n), numparamvariables, numparamfijos, indparamvariables, indparamfijos, valorparamfijos);
end

if numparamvariables == 0 % In the case that there is only one variable to fit, when it is fixed no fitting is needed, only Chi2
    numData=numel(ydata); 
    numParamVariables = 1; % in this case only one parameter is variable in the fit
    ymodel =FUN(valorparamfijos, xdata);
    chi2=sum(((ydata-ymodel)./yerr).^2)/(numData-numParamVariables);
else
    [~, allParam, ~, ~, ~, ~, ~, jacob_mat]...
        = ajusta_lsqnonlin(FUN, guess, LB ,UB ,[], xdata, ydata, yerr, indparamvariables, indparamfijos, valorparamfijos);
    [chi2, ~, ~] = ajusta_computeuncertainties (FUN, xdata, ydata, yerr, allParam, indparamvariables, jacob_mat); %Simplemente calcula el modelo y las incertidumbres del ajuste
end

function chi2 = fitcore_global_SPA (FUN, numParam, dataCell, varGlobals_table, varFix_table, paramlibre, paramsIni_table, valorLB_table, valorUB_table)
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

if isempty(find(varFix == false)) % In the case that there is only one variable to fit, when it is fixed no fitting is needed, only Chi2
    
    chi2 = chi2_global_noFit(dataCell, varGlobals, varFix, paramsIni, valorLB, valorUB, FUN);
    
else % Compute chi2 with fit
    [indparam, fitparam, allParam, indparamvariables, indparamfijos, jacob_mat, dsid]...
        = ajusteGlobalGral_gui(dataCell, varGlobals', paramsIni', varFix', valorLB', valorUB', FUN);

    Gdata=cell2mat(dataCell);
    xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
    ydata=Gdata(:,2);
    yerr=Gdata(:,3);

    [chi2, ~, ~] = ajusta_computeuncertainties_global (FUN, dsid, xdata, ydata, yerr, indparamvariables, indparamfijos, fitparam, allParam, indparam, jacob_mat); %Simplemente calcula el modelo y las incertidumbres del ajuste
end

function chi2 = chi2_global_noFit(dataCell, varGlobals, varFix, paramsIni, valorLB, valorUB, FUN)
% Reestructure the data to compute only Chi2 without fitting, with all the
% variables fixed

    Gdata=cell2mat(dataCell);
    xdata=Gdata(:,1);   % Los datos están estructurados siempre en {x,y,error}
    ydata=Gdata(:,2);
    yerr=Gdata(:,3);
    numCurvas=numel(dataCell);          % Nº de curvas
    numParams=numel(paramsIni);
    numPtosxCurva=zeros(numCurvas,1);   % Nº puntos por curva
    numParamsxCurva=numel(varGlobals);  % Nº de parámetros x curva (a,b,c)
    contadorPtos=0;
    for k=1:numCurvas % Calcula dsid
        dataCurva=dataCell{k};
        numPtosxCurva(k)=size(dataCurva,1);
        desde=contadorPtos+1;
        hasta=contadorPtos+numPtosxCurva(k);
        dsid(desde:hasta,1)=k;
        contadorPtos=hasta;
    end
    indparam = zeros(numCurvas,numParamsxCurva); 
    contadorInd=1; %Contador indices de los parámetros que van a ser ajustados
    contadorFil=1; %Contador de filas de paramsIniFit
    % Calcula indParamsGlobal2Indep, de forma que cada fila contiene los parámetros de ajuste de cada curva
    for aj03=1:numParamsxCurva 
        isGlobal=varGlobals(aj03); 
        switch isGlobal
            case 0
            indparam(:,aj03)=contadorFil:1:contadorFil+numCurvas-1;
            contadorFil=contadorFil+numCurvas;
            case 1
            indparam(:,aj03)=contadorFil*ones(numCurvas,1);
            contadorFil=contadorFil+1;
        end
        contadorInd=contadorInd+numCurvas;
    end
    % Inicializa variables para crear los índices de parámetros fijos y variables
    x = [xdata, dsid];
    numParamVariables=1; % In this case there is only one parameter in the fitting
    numData=numel(ydata);
    FUN_function = str2func(FUN);
    ymodel = FUN_function(paramsIni, indparam, x); % ParamsIni contains the value of the fixed variables
    chi2 = sum(((ydata-ymodel)./yerr).^2)/(numData-numParamVariables);