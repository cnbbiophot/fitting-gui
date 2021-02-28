function PCHnorm=fitfcn_PCHsingleSpecies_OPC1order_sum_global (param, indParam, x)

% PCH for ONE species contained corrected in FIRST order in the article 
% "TPhoton Counting Histogram: One-Photon Excitation". Bo Huang, 
% Thomas D. Perroud, and Richard N. Zare, 2004 that corrects the PCH in 
% "The Photon Counting Histogram in Fluorescence Fluctuation Spectroscopy". 
% Y. Chen et al. 1999
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=averageN
% param(2)=epsilon
% param(3)=F
% param(4)=Q

% agv, 17jun2020

tdata=x(:,1);
dsid=x(:,2); %Indica el nº de curva de cada set de datos

number_sets = length(indParam(:,1));

averageN=param(indParam(:,1));    % Parámetros
epsilon=param(indParam(:,2));
F=param(indParam(:,3));       
Q=param(indParam(:,4));  

k = tdata; % is a vector (1,2..., maxK, 1,2..., maxK, 1,2..., maxK , ...)
length_K = length(k)/number_sets; % Length of photon couns in signals

max_sum = 10000;
p1 = zeros(1,length_K*number_sets);

for ds = 1:number_sets %for the global fit
    % k=0
        k=0;
        i = 2:max_sum;
        sum_term = (-1).^i.*epsilon(ds).^i./(factorial(i).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(length_K*(ds-1)+1) = 1./(2*sqrt(2).*Q(ds).*(1+F(ds))^2)...
            .*(sum(sum_term) - (1+F(ds))*epsilon(ds)) + 1;
    % k=1
        k=1;
        i = 2:max_sum;
        sum_term = (-1).^i.*epsilon(ds).^i./(factorial(i).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(length_K*(ds-1)+2) = -1./(2*sqrt(2)*Q(ds)*(1+F(ds))^2)...
        .*(sum(sum_term) - (1+F(ds))*epsilon(ds));
    % k>=2
    for k = 3:length_K
        real_k = k-1;
        i = real_k:max_sum;
        sum_term = (-1).^i.*epsilon(ds).^i./(factorial(i-real_k).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(length_K*(ds-1)+k) = (-1)^real_k/factorial(real_k)./(2*sqrt(2)*Q(ds)*(1+F(ds))^2)...
        .*sum(sum_term);
    end
end

N = 1E4; % max number of particles to tale into account. Should be higher than 1E3 since there are only 
pN = zeros(N,length_K*number_sets);

for ds = 1:number_sets
    for i = 1:N % We compute all the convolutions in advance
        if i > 2
            convp1pN = conv(p1(length_K*(ds-1)+1:length_K*ds) , pN(i-1,length_K*(ds-1)+1:length_K*ds));
            pN(i,length_K*(ds-1)+1:length_K*ds) = convp1pN(1:length_K);   
        elseif i == 2
            pN(2,length_K*(ds-1)+1:length_K*ds) = p1(length_K*(ds-1)+1:length_K*ds);
        elseif i == 1
            pN(1,length_K*(ds-1)+1) = 1;
        end
    end
end

vectorParticles = 0:N-1;

for ds = 1:number_sets
    PoisDiffusion(ds,:) = exp(-averageN(ds)).*averageN(ds).^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0
end

PCH = zeros(1,length_K*number_sets);
PCHnorm = zeros(length_K*number_sets,1);
for ds = 1:number_sets
        PCH(length_K*(ds-1)+1:length_K*ds) = sum(pN(:,length_K*(ds-1)+1:length_K*ds).*PoisDiffusion(ds,:)')';
        PCHnorm(length_K*(ds-1)+1:length_K*ds,1) = PCH(length_K*(ds-1)+1:length_K*ds)/ sum(PCH(length_K*(ds-1)+1:length_K*ds));
end