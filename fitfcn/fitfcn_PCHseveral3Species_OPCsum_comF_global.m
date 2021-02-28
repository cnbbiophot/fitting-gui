function PCHnorm=fitfcn_PCHseveral3Species_OPCsum_comF_global (param, indParam, x)

% PCH for TWO species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999
%
% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=averageN1
% param(2)=averageN2
% param(3)=averageN3
% param(4)=epsilon1
% param(5)=epsilon2
% param(6)=epsilon3
% param(7)=F
% param(8)=Q

% agv, 15may2020
% agv, 1jun2020

tdata=x(:,1);
dsid=x(:,2); %Indica el nº de curva de cada set de datos

number_sets = length(indParam(:,1));

k = tdata; % is a vector (1,2..., maxK, 1,2..., maxK, 1,2..., maxK , ...)
length_K = length(k)/number_sets; % Length of photon couns in signals

averageN1=param(indParam(:,1));    % Parámetros
averageN2=param(indParam(:,2)); 
averageN3=param(indParam(:,3));
epsilon1=param(indParam(:,4));  
epsilon2=param(indParam(:,5));
epsilon3=param(indParam(:,6));
F=param(indParam(:,7));
Q=param(indParam(:,8));

PCHnorm = simulationPCH_severalSpecies_global([epsilon1 epsilon2 epsilon3],[averageN1 averageN2 averageN3], ...
    F, Q, length_K, number_sets);

end

function PCH_total_norm = simulationPCH_severalSpecies_global(epsilon,averageN,F,Q,length_K, number_sets)

    numberSpecies = length(epsilon(1,:)); % As many species as we input in the function
    
    % Generate every individual probability followong (26) they should be
    % independent so we can convolve them
    PCH_species = zeros(numberSpecies, length_K*number_sets);
    for S = 1:numberSpecies
        PCH_species(S,:) = generatePCH_severalSpecies_global(length_K, epsilon(:,S), averageN(:,S), F, Q, number_sets);        
    end

    PCH_total = PCH_species(1,:); %Convolution to account for several species
    for S = 2:numberSpecies
        for ds = 1:number_sets
            convPCH = conv(PCH_species(S,length_K*(ds-1)+1:length_K*ds),PCH_total(length_K*(ds-1)+1:length_K*ds));
            PCH_total(length_K*(ds-1)+1:length_K*ds) = convPCH(1:length_K);
        end
    end

    PCH_total_norm = zeros(length_K*number_sets,1);
    for ds = 1:number_sets
            PCH_total_norm(length_K*(ds-1)+1:length_K*ds,1) = PCH_total(length_K*(ds-1)+1:length_K*ds)/ sum(PCH_total(length_K*(ds-1)+1:length_K*ds));
    end
end

function PCH_total = generatePCH_severalSpecies_global(length_K, epsilon, averageN, F, Q, number_sets)

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

N = 10000; % max number of particles, it does not matter, since there are only 
% about 1.5 particles in the focus large number of partices dont make the
% difference
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

PCH_total = zeros(1,length_K*number_sets);
for ds = 1:number_sets % Not normalized PCH of one set
        PCH_total(length_K*(ds-1)+1:length_K*ds) = sum(pN(:,length_K*(ds-1)+1:length_K*ds).*PoisDiffusion(ds,:)')';
end

end
