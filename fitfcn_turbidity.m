function y=fitfcn_turbidity (param, lambda)

% y=fitfcn_turbidity (param, lambda)
% Ajuste de la turbidez para restarla de un espectro
%  y=K0*lambda^-k+B
%
% Siguiendo CAST97
%
% lambda es la longitud de onda
%
% B el background
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=K0
% param(2)=k
% param(3)=B
%
% jri 7Oct15

  y = param(1).*lambda.^(-param(2))+param(3);
 
 