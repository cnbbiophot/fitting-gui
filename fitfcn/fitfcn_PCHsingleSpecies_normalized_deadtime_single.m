function PCH_total_nodead=fitfcn_PCHsingleSpecies_normalized_deadtime_single (param, tdata)

% PCH for ONE species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999 with the correction for deadtime equation (7) from "The Photon Counting Histogram in Fluorescence 
% Fluctuation Spectroscopy with Non-Ideal Photodetectors" Hillesheim 2003
% https://www.sciencedirect.com/science/article/pii/S0006349503746220

% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=averageN
% param(2)=epsilon
% param(3)=deadT

% agv, 15may2020

warning('k data of the measurement MUST go from k=2')

k = tdata;
maxK = length(k); % For the moment we dont take into account the k=0 value

averageN=param(1);    % Parámetros
epsilon=param(2);
deadT=param(3);

PCH_total_nodead = zeros(maxK,1); % We set the vector to hold the value k=0

    for kp = 1:maxK
        
        dum_summand1 = generatePCH ( kp, epsilon*(1-kp*deadT), averageN);
        summand1 = sum(dum_summand1);

        dum_summand2 = generatePCH ( kp-1 , epsilon*(1-(kp-1)*deadT), averageN);
        summand2 = sum(dum_summand2);
        
        PCH_total_nodead(kp) = summand1 - summand2;
    end
    
    PCH_total_nodead = PCH_total_nodead / sum(PCH_total_nodead);
    
end

function PCHnorm = generatePCH (maxK, epsilon, averageN)

    if maxK > 0
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
                pN(1,1) = 0; % Probability of having 1 photons with 0 particles (we start with k=1, if we would start with k=0 then it would have to be 1)
            end
            %pN(i,:) = pN(i,:) / trapz(pN(i,:));
        end

        vectorParticles = 0:N-1;
        PoisDiffusion = exp(-averageN).*averageN.^(vectorParticles)./...
                factorial(vectorParticles);
        PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0

        PCH = zeros(1,maxK);
        for k = 1:maxK
            PCH(k) = sum(pN(:,k).*PoisDiffusion');
        end

        PCHnorm = PCH';
    
    else
        PCHnorm = zeros(1,maxK);
        PCHnorm(1) = 1;
    end

end
