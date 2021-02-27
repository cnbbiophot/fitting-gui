 function [yFit, J] = err_ajsGlobalGral(valorParamVarIni, data, numPtosxCurva, varGlobales, indParamsGlobal2Indep, indparamvariables, indparamfijos, valorparamfijos, varFix, funcAjuste);

% Funci�n de ajuste general para  los datos experimentales. Minimiza la diferencia entre los valores experimentales y los te�ricos. 
% Llama a la funci�n general que calcula el jacobiano y a la funci�n que calcula la ecuaci�n a ajustar

% Unai, 18May2015
% agv, 24-Oct-2018

funcAjusteHandle = str2func(funcAjuste); %Handle de la funci�n de ajuste
funcJacobiano=strcat(funcAjuste,'_jac'); %Nombre de la funci�n que calcula el jacobiano de la funci�n de ajuste (funcAjuste) 
funcJacobianoHandle = str2func(funcJacobiano); %Handle de la funci�n del jacobiano

paramsIni(indparamvariables) = valorParamVarIni;
paramsIni(indparamfijos) = valorparamfijos;

paramsIni = paramsIni'; % Hay que hacer esto para que coja los par�metros correctamente

yReal=data(:,3);    % It uses the matrix of the form [x,dsid,yExp,yError]
yErr=data(:,4);
yTeorica=funcAjusteHandle(paramsIni, indParamsGlobal2Indep, data);

yFit=(yTeorica-yReal)./yErr; %Funci�n objetivo, incluye el error en cada punto. La minimizar� lsqnonlin.

numVarFit=numel(varGlobales);
numCurvas=numel(numPtosxCurva);
numVarGlobales=numel(find(varGlobales==1));

dimJacobiano=[size(data,1),numVarGlobales+(numVarFit-numVarGlobales)*numCurvas-numel(indparamfijos)];
J=ajusteGlobalGral_jacobiano(dimJacobiano, paramsIni, indParamsGlobal2Indep,...
    data, varGlobales, numCurvas, numPtosxCurva, varFix, funcJacobianoHandle);

