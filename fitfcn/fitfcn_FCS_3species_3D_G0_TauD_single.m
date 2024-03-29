function G=fitfcn_FCS_3species_3D_G0_TauD_single (param, tdata)

% G=fitfcn_FCS_3DTauD (param, tdata)
% Función de autocorrelación con un modelo de difusion 3D single photon 
%  G = G0*(f1/((1+4*D1*tau/omega0^2)*sqrt(1+4*D1*tau/z0^2))+(1-f1)/((1+4*D2*tau/omega0^2)*sqrt(1+4*D2*tau/z0^2)))
%
% SP = z0/omega0
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=G0
% param(2)=tauD1 [ms]
% param(3)=tauD2 [ms]
% param(4)=tauD3 [ms]
% param(5)=SP
% param(6)=f1
% param(7)=f2
%
% agv 14Dic2018

tau = tdata;

G0 = param(1);
tauD1 = param(2);
tauD2 = param(3);
tauD3 = param(4);
SP = param(5);
f1 = param(6);
f2 = param(7);

 G = G0*(f1./((1+tau/tauD1).*sqrt(1+tau/SP^2/tauD1))+(f2)./((1+tau/tauD2).*sqrt(1+tau/SP^2/tauD2))...
     +(1-f1-f2)./((1+tau/tauD3).*sqrt(1+tau/SP^2/tauD3)));
 