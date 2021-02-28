function fVes=fitfcn_partitionCoefficient_M3_melo11_maxBind_single(param, tdata)
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
% param(6)=mu(fix)
% param(7)=W(fix)
% param(8)=Xi_max
% 
% jri 16May2018
% agv 04Dic2018

% Si no se conoce PTotal hay que trabajar con las fracciones de proteína
% libre (fProt_M3) que se obtienen de los datos de FCS.

% Define the parameters as a function of the input parameters
Kp = param(1);
freeProbe = param(2);
q = param(3);
PTotal = param(4);
VTotal = param(5);
mu = param(6);
W = param(7);
Xi_max = param(8);

LAccessible = tdata;

%Dado Kp, Laccessible y Ptotal obtengo el valor de k_average (número de proteínas asociadas a cada vesícula)
%Eso me permite calcular la probabilidad P de que una vesícula lleve k proteínas

LTotal=2*LAccessible; % Ltotal: Concentración de lípido total 

xm=(Kp.*LAccessible*Xi_max)./(W+Kp.*LAccessible); % Mole fraction of labeled-proteins bound to vesicles
k_average=xm.*PTotal.*mu./LTotal;

NAvogadro=6.022E23;
NVes_total=LTotal*VTotal*NAvogadro/mu;
NProt_total=PTotal*VTotal*NAvogadro;

P_k_k_average=@(x, x_average) (x_average.^x.*exp(-x_average))./factorial(x); % Definición de la probabilidad de Poisson

%P=prob_vesicle_with_k_proteins (k, k_average);
%Calculo la probabilidad de encontrar una vesícula sin proteínas (k=0)
%cuando, en promedio, las vesículas llevan <k> proteínas
P0_k_average=P_k_k_average(0, k_average);
%Nves_prot es el número de vesículas que llevan al menos una proteína
NVes_prot=NVes_total.*(1-P0_k_average);
%Por tanto, Nprot_free es el número de proteínas libres
NProt_free=NProt_total*(1-xm);

k = 1:100;
for i=1:numel(LAccessible) % Each different value of LAccessible has a different 
    sumP(i,1) = sum(k.^2.*P_k_k_average(k, k_average(i)));
end

NDye=freeProbe*NProt_total;
fVes=sumP.*NVes_total./(NDye*q^2+NProt_free+sumP.*NVes_total);