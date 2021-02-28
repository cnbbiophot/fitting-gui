function G=fitfcn_FCS_3D_G0_D_single (param, tdata)

% Función de autocorrelación con un modelo de difusion 3D single photon 
%  G = G0*((1+tdata*4*D/Omega0^2).^-1).*((1+tdata*4*D/z0^2).^-0.5)
%
% SP es el parámetro estructural: S=z0/Omega0
% tdata en s
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=G0
% param(2)=D [um^2/s]
% param(3)=z0 [um]
% param(4)=omega0 [um]

% agv, 12Dic2018

G0=param(1);    % Parámetros
D=param(2);
z0=param(3);
omega0=param(4);

G = G0.*(1+tdata*4.*D./omega0.^2).^(-1).*(1+tdata*4.*D./z0.^2).^(-0.5);
 