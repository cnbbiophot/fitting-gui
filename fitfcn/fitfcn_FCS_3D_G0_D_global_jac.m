function derivParciales=fitfcn_FCS_3D_G0_D_global_jac(params, indParams, data)

% Calcula las derivadas parciales para la función fitfcn_FCS_3DTauD, con respecto a las variables de ajuste (N, tauD, s):
%  G = G0*((1+tdata*4*D/Omega0^2).^-1).*((1+tdata*4*D/z0^2).^-0.5)
% Cada columna corresponde a una variable de ajuste

% agv, 19Oct2018

derivParciales=zeros(size(data,1),3);
G=fitfcn_FCS_3D_G0_D_global (params, indParams, data);
tdata=data(:,1);
dsid=data(:,2); %Indica a qué curva corresponde cada dato
yErr=data(:,4);
paramsCol=params(indParams); %Parámetros de cada curva (por cada columna)
G0=paramsCol(:,1);
G0vector=G0(dsid); %Variable G0 para cada punto del set de datos
D=paramsCol(:,2);
Dvector=D(dsid); %Variable D para cada punto del set de datos
z0=paramsCol(:,3);
z0vector=z0(dsid); %Variable z0 para punto del set de datos
omega0=paramsCol(:,4);
omega0vector=omega0(dsid); %Variable Omega0 para punto del set de datos

% d(G)/d(G0)
derivParciales(:,1)=G./G0vector./yErr;
% d(G)/d(D)
derivParciales(:,2)=-G.*((1+tdata*4.*Dvector./omega0vector.^2).^(-1).*tdata*4./omega0vector.^2+tdata*2./z0vector.^2.*(1+tdata*4.*Dvector./z0vector.^2).^(-1))./yErr; 
% d(G)/d(z0)
derivParciales(:,3)=G.*tdata*4.*Dvector./z0vector.^3.*(1+tdata*4.*Dvector./z0vector.^2).^(-1)./yErr; 
% d(G)/d(Omega0)
derivParciales(:,4)=G.*tdata*8.*Dvector./omega0vector.^3.*(1+tdata*4.*Dvector./omega0vector.^2).^(-1)./yErr; 
