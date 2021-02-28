function PCHnorm=fitfcn_PCH_Boltzmann_dist_sum_comF_single (param, tdata)

% PCH for TWO species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999
%
% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=totalAb [M]
% param(2)=boundAb
% param(3)=x_average
% param(4)=epsilon1
% param(5)=F
% param(6)=omega0 [um]
% param(7)=z0 [um]


% agv, 17feb2021

k = tdata;  % Photon counts from the histogram
length_K = length(k);

concAb=param(1);
boundAb=param(2);
x_average=param(3);
epsilon1=param(4);
F=param(5);
omega0=param(6);
z0=param(7);


volumefocus = (pi/2)^(3/2) * omega0^2 * z0 * 1e-15 * (1+F); % in L
totalAb = concAb*volumefocus*6.022E23;

freeAb = totalAb - boundAb;

averageN1= freeAb + boundAb.*P_k_Boltzmann(1, x_average);    % Parámetros
averageN2= boundAb.*P_k_Boltzmann(2, x_average) ;
averageN3= boundAb.*P_k_Boltzmann(3, x_average) ;
averageN4= boundAb.*P_k_Boltzmann(4, x_average) ;
averageN5= boundAb.*P_k_Boltzmann(5, x_average) ;
averageN6= boundAb.*P_k_Boltzmann(6, x_average) ;
averageN7= boundAb.*P_k_Boltzmann(7, x_average) ;
averageN8= boundAb.*P_k_Boltzmann(8, x_average) ;
averageN9= boundAb.*P_k_Boltzmann(9, x_average) ;

epsilon2=epsilon1*2;
epsilon3=epsilon1*3;
epsilon4=epsilon1*4;
epsilon5=epsilon1*5;
epsilon6=epsilon1*6;
epsilon7=epsilon1*7;
epsilon8=epsilon1*8;
epsilon9=epsilon1*9;

Q = 1; % Fix Q parameter

PCHnorm = simulationPCH_severalSpecies([epsilon1 epsilon2 epsilon3 epsilon4 epsilon5 epsilon6 epsilon7 epsilon8 epsilon9 ],...
    [averageN1 averageN2 averageN3 averageN4 averageN5 averageN6 averageN7 averageN8 averageN9 ], F, Q, length_K);

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

function P_k_Boltzmann = P_k_Boltzmann(x, x_average) 
        P_k_Boltzmann = x_average.^x ./ (x_average + 1).^(x + 1);
end