function G=fitfcn_FCS_2species_3D_G0_D_global (param, indParam, x)

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
% agv 17Dic2018

tdata = x(:,1);
tau = tdata;
dsid = x(:,2); %Indica el nº de curva de cada set de datos

G0 = param(indParam(:,1));
D1 = param(indParam(:,2));
D2 = param(indParam(:,3));
omega0 = param(indParam(:,4));
z0 = param(indParam(:,5));
f1 = param(indParam(:,6));


 G = G0(dsid).*(f1(dsid)./((1+4*D1(dsid).*tau./omega0(dsid).^2).*sqrt(1+4*D1(dsid).*tau./z0(dsid).^2))...
     +(1-f1(dsid))./((1+4*D2(dsid).*tau./omega0(dsid).^2).*sqrt(1+4*D2(dsid).*tau./z0(dsid).^2)));
 