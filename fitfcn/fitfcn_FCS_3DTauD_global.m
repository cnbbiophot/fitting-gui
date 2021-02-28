function G=fitfcn_FCS_3DTauD_global (param, indParam, x)
           
% G=fitfcn_FCS_3DTauD (param, tdata)
% Funci�n de autocorrelaci�n con un modelo de difusion 3D single photon 
%  G = (1/N)*((1+tdata/tauD).^-1).*((1+tdata/(tauD*SP^2)).^-0.5)
%
% SP es el par�metro estructural: S=z0/w0
% tdata en s!!
% tauD en ms!!
%
% Definici�n de los par�metros actual: param(n)=Parameter Units
% Definici�n de los par�metros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=N
% param(2)=tauD [ms]
% param(3)=SP
%
% jri 4Feb15
% Unai 6Mar15, modificado para ajuste global

tdata=x(:,1);
dsid=x(:,2); %Indica el n� de curva de cada set de datos
N=param(indParam(:,1));
tauD=param(indParam(:,2));
s=param(indParam(:,3));

G = ((1+tdata./(tauD(dsid)*1E-3)).^-1).*((1+tdata./(tauD(dsid).*1E-3.*s(dsid).^2)).^-0.5)./N(dsid);
 