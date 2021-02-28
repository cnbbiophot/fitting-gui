function PCHnorm=fitfcn_PCHseveral2Species_normalized_global (param, indParam, x)

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
% param(3)=epsilon1
% param(4)=epsilon2

% agv, 15may2020
% agv, 1jun2020

tdata=x(:,1);
dsid=x(:,2); %Indica el nº de curva de cada set de datos

number_sets = length(indParam(:,1));

k = tdata; % is a vector (1,2..., maxK, 1,2..., maxK, 1,2..., maxK , ...)
maxK = length(k)/number_sets; % For the moment we dont take into account the k=0 value

averageN1=param(indParam(:,1));    % Parámetros
averageN2=param(indParam(:,2)); 
epsilon1=param(indParam(:,3));  
epsilon2=param(indParam(:,4));

PCHnorm = simulationPCH_severalSpecies_global([epsilon1 epsilon2],[averageN1 averageN2], maxK, number_sets, k);

end

function PCH_total_norm = simulationPCH_severalSpecies_global(epsilon,averageN,maxK, number_sets, k)

    numberSpecies = length(epsilon(1,:)); % As many species as we input in the function
    
    % Generate every individual probability followong (26) they should be
    % independent so we can convolve them
    PCH_species = zeros(numberSpecies, maxK*number_sets);
    for S = 1:numberSpecies
        PCH_species(S,:) = generatePCH_severalSpecies_global(maxK, epsilon(:,S), averageN(:,S), number_sets, k);        
    end

    PCH_total = PCH_species(1,:); %Convolution to account for several species
    for S = 2:numberSpecies
        for ds = 1:number_sets
            convPCH = conv(PCH_species(S,maxK*(ds-1)+1:maxK*ds),PCH_total(maxK*(ds-1)+1:maxK*ds));
            PCH_total(maxK*(ds-1)+1:maxK*ds) = convPCH(1:maxK);
        end
    end

    PCH_total_norm = zeros(maxK*number_sets,1);
    for ds = 1:number_sets
            PCH_total_norm(maxK*(ds-1)+1:maxK*ds,1) = PCH_total(maxK*(ds-1)+1:maxK*ds)/ sum(PCH_total(maxK*(ds-1)+1:maxK*ds));
    end
end

function PCH_total = generatePCH_severalSpecies_global(maxK, epsilon, averageN, number_sets, k)

% Compute 1 particle probability with formula (16)
p1 = zeros(1,maxK*number_sets);

% we elliminate the for to introduce the dsid 
for ds = 1:number_sets %for the global fit
    PoisInt = @(x_fun) gammainc(abs(epsilon(ds)).*exp(-2*x_fun.^2),k(maxK*(ds-1)+1:maxK*ds)) .* gamma(k(maxK*(ds-1)+1:maxK*ds)); % Careful with the order of the input factors
    p1(maxK*(ds-1)+1:maxK*ds) = 1./(factorial(k(maxK*(ds-1)+1:maxK*ds))) .* integral(PoisInt,0,Inf,'ArrayValued',true); 
end

N = 10000; % max number of particles, it does not matter, since there are only 
% about 1.5 particles in the focus large number of partices dont make the
% difference
pN = zeros(N,maxK*number_sets);

for ds = 1:number_sets
    for i = 1:N % We compute all the convolutions in advance
        if i > 2
            convp1pN = conv(p1(maxK*(ds-1)+1:maxK*ds) , pN(i-1,maxK*(ds-1)+1:maxK*ds));
            pN(i,maxK*(ds-1)+1:maxK*ds) = convp1pN(1:maxK);   
        elseif i == 2
            pN(2,maxK*(ds-1)+1:maxK*ds) = p1(maxK*(ds-1)+1:maxK*ds);
        elseif i == 1
            pN(1,maxK*(ds-1)+1) = 0; % Probability of having 1 photons with 0 particles (we start with k=1, if we would start with k=0 then it would have to be 1)
        end
    end
end

vectorParticles = 0:N-1;

for ds = 1:number_sets
    PoisDiffusion(ds,:) = exp(-averageN(ds)).*averageN(ds).^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0
end

PCH_total = zeros(1,maxK*number_sets);
for ds = 1:number_sets % Not normalized PCH of one set
        PCH_total(maxK*(ds-1)+1:maxK*ds) = sum(pN(:,maxK*(ds-1)+1:maxK*ds).*PoisDiffusion(ds,:)')';
end

end
