function PCHnorm=fitfcn_PCHseveral2Species_normalized_single (param, tdata)

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

% agv, 29may2020

k = tdata;
maxK = length(k); % For the moment we dont take into account the k=0 value

averageN1=param(1);    % Parámetros
averageN2=param(2);
epsilon1=param(3);
epsilon2=param(4);

PCHnorm = simulationPCH_severalSpecies([epsilon1 epsilon2],[averageN1 averageN2], maxK);

end

function PCH_total_norm = simulationPCH_severalSpecies(epsilon,averageN,maxK)

    numberSpecies = length(epsilon); % As many species as we input in the function
    
    % Generate every individual probability followong (26) they should be
    % independent so we can convolve them
    PCH_species = zeros(numberSpecies, maxK);
    for S = 1:numberSpecies
        PCH_species(S,:) = generatePCH(maxK, epsilon(S), averageN(S));        
    end

    PCH_total = PCH_species(1,:);
    for S = 2:numberSpecies
        convPCH = conv(PCH_species(S,:),PCH_total);
        PCH_total = convPCH(1:maxK);
    end

    PCH_total_norm = PCH_total' / sum(PCH_total);
end

function PCH = generatePCH(maxK, epsilon, averageN)

    % Compute 1 particle probability with formula (16)
    
    p1 = zeros(1,maxK);
    for k = 1:maxK
        PoisInt = @(x) gammainc(epsilon.*exp(-2*x.^2),k) * gamma(k); % Careful with the order of the input factors
        p1(k) = 1/(factorial(k)) .* integral(PoisInt,0,Inf); 
    end

    N = 10000; % max number of particles, it does not matter, since there are only 
    % about 1.5 particles in the focus large number of partices dont make the
    % difference
    pN = zeros(N,maxK);

    for i = 1:N % We compute all the convolutions in advance
        if i > 2
            convp1pN = conv(p1,pN(i-1,:));
            pN(i,:) = convp1pN(1:maxK);   
        elseif i == 2
            pN(2,:) = p1; % First convolution
        elseif i == 1
            pN(1,1) = 0; % Probability of having 1 photons with 0 particles (we start with k=1, if we woudl start with k=0 then it would have to be 1)
        end
    end
    
    vectorParticles = 0:N-1;
    PoisDiffusion = exp(-averageN).*averageN.^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0

    PCH = zeros(1,maxK);
    for k = 1:maxK
        PCH(k) = sum(pN(:,k).*PoisDiffusion');
    end
end