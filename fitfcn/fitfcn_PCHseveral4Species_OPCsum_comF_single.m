function PCHnorm=fitfcn_PCHseveral4Species_OPCsum_comF_single (param, tdata)

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
% param(4)=averageN4
% param(5)=epsilon1
% param(6)=epsilon2
% param(7)=epsilon3
% param(8)=epsilon4
% param(9)=F
% param(10)=Q

% agv, 17jun2020

k = tdata;  % Photon counts from the histogram
length_K = length(k);

averageN1=param(1);    % Parámetros
averageN2=param(2);
averageN3=param(3);
averageN4=param(4);
epsilon1=param(5);
epsilon2=param(6);
epsilon3=param(7);
epsilon4=param(8);
F=param(9);
Q=param(10);

PCHnorm = simulationPCH_severalSpecies([epsilon1 epsilon2 epsilon3 epsilon4],[averageN1 averageN2 averageN3 averageN4],...
    F, Q, length_K);

end

function PCH_total_norm = simulationPCH_severalSpecies(epsilon,averageN, F, Q,length_K)

    numberSpecies = length(epsilon); % As many species as we input in the function
    
    % Generate every individual probability followong (26) they should be
    % independent so we can convolve them
    PCH_species = zeros(numberSpecies, length_K);
    for S = 1:numberSpecies
        PCH_species(S,:) = generatePCH(length_K, epsilon(S), averageN(S), F, Q);        
    end

    PCH_total = PCH_species(1,:);
    for S = 2:numberSpecies
        convPCH = conv(PCH_species(S,:),PCH_total);
        PCH_total = convPCH(1:length_K);
    end

    PCH_total_norm = PCH_total' / sum(PCH_total);
end

function PCH = generatePCH(length_K, epsilon, averageN, F, Q)

    % Compute 1 particle probability with formula (16)
    
    max_sum = 10000;
    p1 = zeros(1,length_K);
    % k=0
        k=0;
        i = 2:max_sum;
        sum_term = (-1).^i.*epsilon.^i./(factorial(i).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(1) = 1./(2*sqrt(2)*Q*(1+F)^2)...
            .*(sum(sum_term) - (1+F)*epsilon)...
            + 1;
    % k=1
        k=1;
        i = 2:max_sum;
        sum_term = (-1).^i.*epsilon.^i./(factorial(i).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(2) = -1./(2*sqrt(2)*Q*(1+F)^2)...
        .*(sum(sum_term) - (1+F)*epsilon);
    % k>=2
    for k = 3:length_K
        real_k = k-1;
        i = real_k:max_sum;
        sum_term = (-1).^i.*epsilon.^i./(factorial(i-real_k).*sqrt(i).*i);
        sum_term(isnan(sum_term))=0;
        p1(k) = (-1)^real_k/factorial(real_k)./(2*sqrt(2)*Q*(1+F)^2)...
        .*sum(sum_term);
    end

    N = 10000; % max number of particles, it does not matter, since there are only 
    % about 1.5 particles in the focus large number of partices dont make the
    % difference
    pN = zeros(N,length_K);

    for i = 1:N % We compute all the convolutions in advance
        if i > 2
            convp1pN = conv(p1,pN(i-1,:));
            pN(i,:) = convp1pN(1:length_K);   
        elseif i == 2
            pN(2,:) = p1; % First convolution
        elseif i == 1
            pN(1,1) = 1;
        end
    end
    
    vectorParticles = 0:N-1;
    PoisDiffusion = exp(-averageN).*averageN.^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0

    PCH = zeros(1,length_K);
    for k = 1:length_K
        PCH(k) = sum(pN(:,k).*PoisDiffusion');
    end
end