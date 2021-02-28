function G=fitfcn_FCS_gauss3DG0 (param, tdata)

% G=fitfcn_FCS_gauss3DG0 (param, tdata);
% Funci�n de autocorrelaci�n con un modelo de difusion 3D single photon y haz gaussiano
%
% G = G0./((1+4*D*tdata/w0^2).*(1+4*D*tdata/z0^2).^-0.5);
%
% Definici�n de los par�metros actual: param(n)=Parameter Units
% Definici�n de los par�metros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=omega0 [mum]
% param(2)=z0 [mum]
% param(3)=G0 
% param(4)=D [mum^2/s]
%
% jri & GdlH (12nov10)


 G = param(3)./((1+4*param(4)*tdata/param(1)^2).*(1+4*param(4)*tdata/param(2)^2).^-0.5);
 
