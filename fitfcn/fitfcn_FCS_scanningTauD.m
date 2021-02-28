function G=fitfcn_FCS_scanningTauD (param, tdata)

% G=fitfcn_FCS_scanningTauD (param, tdata)
% Funci�n de autocorrelaci�n con un modelo de difusion 2D single photon 
%  G = (1/N)*((1+tdata/tauD).^-.5).*((1+tdata/(tauD*SP^2)).^-0.5)
%
% tdata en s!!
% tauD en ms!!
% SP es el par�metro estructural: S=z0/w0
%
% Definici�n de los par�metros actual: param(n)=Parameter Units
% Definici�n de los par�metros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=N
% param(2)=tauD [ms]
% param(3)=SP
%
% jri 4Feb15

 
 G = ((1+tdata/(param(2)*1E-3)).^-.5).*((1+tdata/(param(2)*1E-3*param(3)^2)).^-0.5)/param(1);
 
 