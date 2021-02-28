function G=fitfcn_FCS_2D_G0_D_global (param, indParam, x)

% Función de autocorrelación con un modelo de difusion 2D single photon en
% GUVs (escaneo de la membrana)
%
%  G = G0*(((1+tdata*4*D/Omega0^2)*(1+tdata*4*D/z0^2)).^-0.5)
%
% tdata en s
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=G0
% param(2)=D [um^2/s]
% param(3)=omega0 [um]
% param(4)=z0 [um]

% agv, 19Oct2018


tdata=x(:,1);   % Parámetro temporal de la curva, datos temporales de la ACF
dsid=x(:,2);    % Indica el nº de curva de cada set de datos
G0=param(indParam(:,1));    % Parámetros
D=param(indParam(:,2));
omega0=param(indParam(:,3));
z0=param(indParam(:,4));

G = G0(dsid).*(((1+tdata*4.*D(dsid)./omega0(dsid).^2).*(1+tdata*4.*D(dsid)./z0(dsid).^2)).^(-0.5));