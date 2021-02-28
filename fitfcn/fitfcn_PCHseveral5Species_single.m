function PCHnorm=fitfcn_PCHseveral5Species_single (param, tdata)

% PCH for ONE species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999
%
% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=conc1 [M]
% param(2)=conc2 [M]
% param(3)=conc3 [M]
% param(4)=conc4 [M]
% param(5)=conc5 [M]
% param(6)=epsilon1
% param(7)=epsilon2
% param(8)=epsilon2
% param(9)=epsilon4
% param(10)=epsilon5
% param(11)=z0 [m]
% param(12)=omega0 [m]

% agv, 29may2020

k = tdata;
maxK = length(k); % For the moment we dont take into account the k=0 value

conc1=param(1);    % Parámetros
conc2=param(2);
conc3=param(3);
conc4=param(4);
conc5=param(5);
epsilon1=param(6);
epsilon2=param(7);
epsilon3=param(8);
epsilon4=param(9);
epsilon5=param(10);
z0=param(11);
omega0=param(12);

PCHnorm = simulationPCH_severalSpecies(omega0,z0,[epsilon1 epsilon2 epsilon3 epsilon4 epsilon5],[conc1 conc2 conc3 conc4 conc5], maxK);

end

function PCH_total_norm = simulationPCH_severalSpecies(omega0,z0,epsilon,conc,maxK)

    numberSpecies = length(epsilon); % As many species as we input in the function
    
    % Generate every individual probability followong (26) they should be
    % independent so we can convolve them
    PCH_species = zeros(numberSpecies, maxK);
    for S = 1:numberSpecies
        PCH_species(S,:) = generatePCH(maxK, omega0, z0, epsilon(S), conc(S));        
    end

    PCH_total = PCH_species(1,:);
    for S = 2:numberSpecies
        convPCH = conv(PCH_species(S,:),PCH_total);
        PCH_total = convPCH(1:maxK);
    end

    PCH_total_norm = PCH_total' / sum(PCH_total);
end

function PCH = generatePCH(maxK, omega0, z0, epsilon, conc)

    VPSF = 4/3*pi*omega0^2*z0; % Volume of PSF in m^3
    Na = 6.022E23;
    % Compute 1 particle probability with formula (16)
    
    p1 = zeros(1,maxK);
    for k = 1:maxK
        PoisInt = @(x) gammainc(epsilon.*exp(-2*x.^2),k) * gamma(k); % Careful with the order of the input factors
        p1(k) = pi*omega0^2*z0/(VPSF * factorial(k)) .* integral(PoisInt,0,Inf); 
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
    
    averageN = conc * VPSF * 1E3 * Na; % Average number of particles
    vectorParticles = 0:N-1;
    PoisDiffusion = exp(-averageN).*averageN.^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0

    PCH = zeros(1,maxK);
    for k = 1:maxK
        PCH(k) = sum(pN(:,k).*PoisDiffusion');
    end
end