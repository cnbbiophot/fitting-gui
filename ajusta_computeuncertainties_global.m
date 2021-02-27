function [Chi2, alldeltaParam, ymodel] = ajusta_computeuncertainties_global (FUN, dsid, xdata, ydata, yerr, indparamvariables, indparamfijos, fitparam, allparam, indparam, jacob_mat)

% [Chi2 deltaParam ymodel] = ajusta_computeuncertainties (FUN, xdata, ydata, yerr, allparam, indparamvariables, jacob_mat);
% Calcula el Chi2 y las incertidumbres de los par�metros ajustados por ajusta_lsqnonlin

%  FUN es un handle a la funci�n de ajuste
%  paramfit son los resultados del ajuste
%  jacob_mat es una matriz que contiene el jacobiano calculado por Matlab (corresponde a las derivadas de G respecto de los par�metros ajustados)
%
%  chi2 es chi2
%  deltaParam son las incertidumbres del ajusta
%  ymodel es el modelo (s�lo los datos y)
%
% jri 3Feb2015
% agv 02-Nov-2018

FUN_function = str2func(FUN);

jacob_mat=full(jacob_mat);

numParamVariables=numel(fitparam);
alfa=zeros(numParamVariables);
numData=numel(ydata); 

alfa = jacob_mat'*jacob_mat;

C=sqrt(inv(alfa));
for n=1:numParamVariables
    deltaParam(n)=C(n,n);
end

x = [xdata, dsid];
ymodel = FUN_function(allparam', indparam, x);

% El chi2 reducido se calcula dividiendo por el # de grados de libertad 
% (calculado como la diferencia entre el # total de puntos de la correlaci�n 
% y el # de par�metros que queremos calcular)

Chi2=sum(((ydata-ymodel)./yerr).^2);
% Compute Chi2 reduced = Chi2/v
% v = N - m; 
% N (data points), m (parameters)
v = (numData-numParamVariables);
Chi2 = Chi2/v;

numtotalparams = numel(indparamvariables)+numel(indparamfijos);
alldeltaParam = zeros (1, numtotalparams);
alldeltaParam (indparamvariables) = deltaParam;





