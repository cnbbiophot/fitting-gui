function [indParamsGlobal2Indep, paramfitY, allparam, indparamvariables, indparamfijos, J, dsid]...
    =ajusteGlobalGral_gui(dataCell, varGlobals, paramsIni, varFix, valorLB, valorUB, funcAjuste)

% Ajuste de una función con variables globales y con elección de parámetros fijos.
% INPUTS:
% dataCell - Celda COLUMNA en la que cada elemento incluye un set de datos completo. 
% varGlobales - Matriz lógica, que indica qué variables de ajuste son globales. CUIDADO CON EL ORDEN!!! Debe ir acorde con los parámetros de la función de ajuste.
% paramsIni - Guesses iniciales del ajuste, en una columna. CUIDADO CON EL ORDEN Y EL TAMAÑO!!!
% varFix - Vector lógica, que indica qué variables de ajuste son FIJAS (==length(paramsIni)).
% paramsFix - Vector con los parámetros fijos.
% funcAjuste - Nombre de la función que será ajustada con lsqnonlin. La función del jacobiano se llamará igual que la de ajuste, añadiendo '_jac' al final.
% 
% OUTPUTS:
% paramfitY - Parámetros ajustados con lsqnonlin. Misma dimensión que paramsIni. 
% allparams - Vector con todos los parámetros, los fijos y los ajustados, con el mismo orden que paramsIni.
% J - Jacobiano de la función.
%
% EJEMPLO DE USO:
% (Las variables son las extraidas por el programa gui_anaFCS.m)
% load('D:\Usuarios\Arturo\180307 FCS Volume Calibration\B&H data\analisis int10freq1000sin18\180307 633_10nm+t5_m2.mat')
% dataCell{1,1}=Gmean ;
% load('D:\Usuarios\Arturo\180307 FCS Volume Calibration\B&H data\analisis int10freq1000sin18\180307 633_20nm+t5_m1.mat')
% dataCell{2,1}=Gmean;
% load('D:\Usuarios\Arturo\180307 FCS Volume Calibration\B&H data\analisis int10freq1000sin18\180307 633_30nm+t5_m1.mat')
% dataCell{3,1}=Gmean;
% paramsIni = [0.01 0.3 303 0.2 0.2 1]'
% varFix = [0 0 1 0 0 0]'
% varGlobales = [0 1 1 0]'
% funcAjuste = 'fitfcn_FCS_3D_G0_D_z0_omega0'
% [paramfitY, allparam, J]=ajusteGlobalGral(dataCell, varGlobales, paramsIni, varFix, funcAjuste);

% Este programa usa las subfunciones: check_fijoovariable_global,
% ajusta_lsqnonlin_global, err_ajsGlobalGral, ajusteGlobalGral_jacobiano y las dos funciones
% definiendo la función de ajuste y el jacobiano (e.g.: fitfcn_FCS_3D_G0_D_z0_omega0 
% y fitfcn_FCS_3D_G0_D_z0_omega0_jac)

% Unai, 19Feb2015
% agv, 02-Nov-2018

numCurvas=numel(dataCell);          % Nº de curvas
numParams=numel(paramsIni);
numParamsxCurva=numel(varGlobals);  % Nº de parámetros x curva (a,b,c)
dataMat=cell2mat(dataCell);         % Datos a ajustar
numPtos=size(dataMat,1);            % Número de puntos
dsid=zeros(numPtos,1); % Indica el nº de curva de cada set de datos de dataMat
numPtosxCurva=zeros(numCurvas,1);   % Nº puntos por curva
contadorPtos=0;

% Comprobación de errores
if numel(varFix) ~= numel(paramsIni)
    error ('¡varFix y paramsIni tienen que tener la misma longitud!'); 
end
if numel(find(varGlobals==0))*numCurvas + numel(find(varGlobals==1)) ~= numel(paramsIni)
    error('Las variables globales e individuales no se corresponden con los valores iniciales introducidos');
end

for aj01=1:numCurvas % Calcula dsid
    dataCurva=dataCell{aj01};
    numPtosxCurva(aj01)=size(dataCurva,1);
    desde=contadorPtos+1;
    hasta=contadorPtos+numPtosxCurva(aj01);
    dsid(desde:hasta,1)=aj01;
    contadorPtos=hasta;
end

data_all = [dataMat(:,1), dsid, dataMat(:,2), dataMat(:,3)]; % [Xdata, dsid, Ydata, Yerr]
%indParamsGlobal2Indep = Indices de conversión de los parámetros de ajuste.Se le pasa a lsqnonlin para acelerar la conversión del formato global al "independiente".
indParamsGlobal2Indep = zeros(numCurvas,numParamsxCurva); 
contadorInd=1; %Contador indices de los parámetros que van a ser ajustados
contadorFil=1; %Contador de filas de paramsIniFit

% Calcula indParamsGlobal2Indep, de forma que cada fila contiene los parámetros de ajuste de cada curva
for aj03=1:numParamsxCurva 
    isGlobal=varGlobals(aj03); 
    switch isGlobal
        case 0
        indParamsGlobal2Indep(:,aj03)=contadorFil:1:contadorFil+numCurvas-1;
        contadorFil=contadorFil+numCurvas;
        case 1
        indParamsGlobal2Indep(:,aj03)=contadorFil*ones(numCurvas,1);
        contadorFil=contadorFil+1;
    end
    contadorInd=contadorInd+numCurvas;
end

% Inicializa variables para crear los índices de parámetros fijos y variables
indparamvariables=[];
indparamfijos=[];
valorparamfijos=[];
numparamfijos=0;
numparamvariables=0;
guess=[];
LB = [];
UB = [];

% Calcula el índice de parámetros variables y fijos
for n=1:numParams
    
    [LB, UB, numparamvariables, numparamfijos, indparamvariables, indparamfijos, valorparamfijos]=...
        check_fijoovariable_global (varFix(n), paramsIni(n), n, paramsIni, LB, UB, valorLB(n), ...
        valorUB(n), numparamvariables, numparamfijos, indparamvariables, indparamfijos, valorparamfijos);
end

valorParamVarIni = paramsIni(indparamvariables);
JacobianOption = 1;
[paramfitY, allparam, indparamvariables,~,~,~,~,~,J] = ...
    ajusta_lsqnonlin_global_gui(valorParamVarIni, LB, UB, ...
    data_all, indparamvariables, indparamfijos, valorparamfijos,...
    numPtosxCurva, varGlobals, indParamsGlobal2Indep, varFix, funcAjuste);