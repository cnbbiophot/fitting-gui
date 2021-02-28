function fVes=fitfcn_partitionCoefficient_M3_melo11_global(param, indParam, x)
%Ecuación 20 de Melo11
% Kp, freeProbe, W, LAccessible, PTotal, VTotal, mu, q
% fVes_M3 = sumP.*NVes_total./((freeProbe*NProt_total)*q^2+NProt_free+sumP.*NVes_total);
% sumP is a function of Kp
%
% Parameters to fit are: Kp, freeProbe and q
%
% Kp: partition constant
% freeProbe: ratio of free proble in solution
% Laccessible: Concentración de lípido accesible Laccessible=Ltotal/2
% LAccessible == tdata
% Ptotal: Concentración total de ab
% Vtotal Volumen total: 150 uL
% Nprot_total= Número total de proteínas
% Nves_total Número total de vesículas
% mu: número de lípidos por vesícula
% q: brillo de la sonda libre con respecto a la sonda ligada a la proteína
% <k> número de proteínas en promedio asociadas a cada vesícula
%
% param(1)=Kp
% param(2)=freeProbe
% param(3)=q
% param(4)=Ptotal(fix)
% param(5)=Vtotal(fix)
% param(6)=Nprot_total(fix)
% param(7)=mu(fix)
% param(8)=W(fix)
% 
% jri 16May2018
% agv 17Dic2018

% Si no se conoce PTotal hay que trabajar con las fracciones de proteína
% libre (fProt_M3) que se obtienen de los datos de FCS.

tdata=x(:,1);
dsid=x(:,2); %Indica el nº de curva de cada set de datos

% Define the parameters as a function of the input parameters
Kp = param(indParam(:,1));
freeProbe = param(indParam(:,2));
q = param(indParam(:,3));
PTotal = param(indParam(:,4));
VTotal = param(indParam(:,5));
NProt_total = param(indParam(:,6));
mu = param(indParam(:,7));
W = param(indParam(:,8));

LAccessible = tdata;

%Dado Kp, Laccessible y Ptotal obtengo el valor de k_average (número de proteínas asociadas a cada vesícula)
%Eso me permite calcular la probabilidad P de que una vesícula lleve k proteínas

LTotal=2*LAccessible; % Ltotal: Concentración de lípido total 

xm=(Kp(dsid).*LAccessible)./(W(dsid)+Kp(dsid).*LAccessible); % Mole fraction of labeled-proteins bound to vesicles
k_average=xm.*PTotal(dsid).*mu(dsid)./LTotal;

NAvogadro=6.022E23;
NVes_total=LTotal*VTotal(dsid).*NAvogadro/mu(dsid);
NProt_total=PTotal(dsid)*VTotal(dsid).*NAvogadro;

P_k_k_average=@(x, x_average) (x_average.^x.*exp(-x_average))./factorial(x); % Definición de la probabilidad de Poisson

%P=prob_vesicle_with_k_proteins (k, k_average);
%Calculo la probabilidad de encontrar una vesícula sin proteínas (k=0)
%cuando, en promedio, las vesículas llevan <k> proteínas
P0_k_average=P_k_k_average(0, k_average);
%Nves_prot es el número de vesículas que llevan al menos una proteína
NVes_prot=NVes_total.*(1-P0_k_average);
%Por tanto, Nprot_free es el número de proteínas libres
NProt_free=NProt_total(dsid).*(1-xm);

k = 1:100;
for i=1:numel(LAccessible) % Each different value of LAccessible has a different 
    sumP(i,1) = sum(k.^2.*P_k_k_average(k, k_average(i)));
end

NDye=freeProbe(dsid).*NProt_total(dsid);
fVes=sumP.*NVes_total./(NDye(dsid).*q(dsid).^2+NProt_free+sumP.*NVes_total);