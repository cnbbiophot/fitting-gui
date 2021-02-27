 function [yFit, J] = err_ajsGlobalGral(valorParamVarIni, data, numPtosxCurva, varGlobales, indParamsGlobal2Indep, indparamvariables, indparamfijos, valorparamfijos, varFix, funcAjuste);

% Función de ajuste general para  los datos experimentales. Minimiza la diferencia entre los valores experimentales y los teóricos. 
% Llama a la función general que calcula el jacobiano y a la función que calcula la ecuación a ajustar

% Unai, 18May2015
% agv, 24-Oct-2018

funcAjusteHandle = str2func(funcAjuste); %Handle de la función de ajuste
funcJacobiano=strcat(funcAjuste,'_jac'); %Nombre de la función que calcula el jacobiano de la función de ajuste (funcAjuste) 
funcJacobianoHandle = str2func(funcJacobiano); %Handle de la función del jacobiano

paramsIni(indparamvariables) = valorParamVarIni;
paramsIni(indparamfijos) = valorparamfijos;

paramsIni = paramsIni'; % Hay que hacer esto para que coja los parámetros correctamente

yReal=data(:,3);    % It uses the matrix of the form [x,dsid,yExp,yError]
yErr=data(:,4);
yTeorica=funcAjusteHandle(paramsIni, indParamsGlobal2Indep, data);

yFit=(yTeorica-yReal)./yErr; %Función objetivo, incluye el error en cada punto. La minimizará lsqnonlin.

numVarFit=numel(varGlobales);
numCurvas=numel(numPtosxCurva);
numVarGlobales=numel(find(varGlobales==1));

dimJacobiano=[size(data,1),numVarGlobales+(numVarFit-numVarGlobales)*numCurvas-numel(indparamfijos)];
J=ajusteGlobalGral_jacobiano(dimJacobiano, paramsIni, indParamsGlobal2Indep,...
    data, varGlobales, numCurvas, numPtosxCurva, varFix, funcJacobianoHandle);

