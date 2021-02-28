function G=fitfcn_FCS_3DTauD_triplet (param, tdata)

% G=fitfcn_FCS_3DTau_triplet (param, tdata)
% Funci�n de autocorrelaci�n con un modelo de difusion 3D single photon 
%  G =
%  (1/N)*(1-F+F*exp(-tdata/tauTriplet))/(1-F)*((1+tdata/tauD).^-1).*((1+tdata/(tauD*SP^2)).^-0.5)
%
% SP es el par�metro estructural: S=z0/w0
% tdata en s!!
% tauTriplet en ms!!
% tauD en ms!!
%
% Definici�n de los par�metros actual: param(n)=Parameter Units
% Definici�n de los par�metros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=N
% param(2)=F
% param(3)=tauTriplet [ms]
% param(4)=tauD [ms]
% param(5)=SP
%
% jri 4Feb15


 G = (1-param(2)+param(2)*exp(-tdata/(param(3)*1E-3))).*((1+tdata/(param(4)*1E-3)).^-1).*((1+tdata/(param(4)*1E-3*param(5)^2)).^-0.5)/(param(1)*(1-param(2)));
 
