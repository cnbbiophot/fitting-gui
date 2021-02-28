function G=fitfcn_FCS_2species_3D_G0_D_single (param, tdata)

% Función de autocorrelación con un modelo de difusion 3D single photon 
%  G = G0*(f1/((1+4*D1*tau/omega0^2)*sqrt(1+4*D1*tau/z0^2))+(1-f1)/((1+4*D2*tau/omega0^2)*sqrt(1+4*D2*tau/z0^2)))
%
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=G0
% param(2)=D1 [um^2/s]
% param(3)=D2 [um^2/s]
% param(4)=omega0 [um]
% param(5)=z0 [um]
% param(6)=f1
%
% agv 14Dic2018

tau = tdata;

G0 = param(1);
D1 = param(2);
D2 = param(3);
omega0 = param(4);
z0 = param(5);
f1 = param(6);


 G = G0*(f1./((1+4*D1*tau/omega0^2).*sqrt(1+4*D1*tau/z0^2))+(1-f1)./((1+4*D2*tau/omega0^2).*sqrt(1+4*D2*tau/z0^2)));
 