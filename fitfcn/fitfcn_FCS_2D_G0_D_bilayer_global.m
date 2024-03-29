function G=fitfcn_FCS_2D_G0_D_bilayer_global (param, indParam, x)

% Funci�n de autocorrelaci�n con un modelo de difusion 3D single photon 
%  G = G0*((1+tdata*4*D/Omega0^2).^-1)
%
% tdata en s
%
% Definici�n de los par�metros actual: param(n)=Parameter Units
% Definici�n de los par�metros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=G0
% param(2)=D [um^2/s]
% param(3)=omega0 [um]

% agv, 19Oct2018


tdata=x(:,1);   % Par�metro temporal de la curva, datos temporales de la ACF
dsid=x(:,2);    % Indica el n� de curva de cada set de datos
G0=param(indParam(:,1));    % Par�metros
D=param(indParam(:,2));
omega0=param(indParam(:,3));

G = G0(dsid).*(1+tdata*4.*D(dsid)./omega0(dsid).^2).^(-1);