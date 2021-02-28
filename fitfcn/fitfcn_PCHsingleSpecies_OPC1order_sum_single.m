function PCHnorm=fitfcn_PCHsingleSpecies_OPC1order_sum_single (param, tdata)

% PCH for ONE species contained corrected in FIRST order in the article 
% "TPhoton Counting Histogram: One-Photon Excitation". Bo Huang, 
% Thomas D. Perroud, and Richard N. Zare, 2004 that corrects the PCH in 
% "The Photon Counting Histogram in Fluorescence Fluctuation Spectroscopy". 
% Y. Chen et al. 1999
%
% from equations (44-45) in "Fluorescence Correlation Spectroscopy (FCS)
% Technical Manual Appendix II PCH"
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=averageN
% param(2)=epsilon
% param(3)=F
% param(4)=Q

% agv, 13jun2020

k = tdata;  % Photon counts from the histogram
length_K = length(k);

averageN=param(1);
epsilon=param(2);
F=param(3);       
Q=param(4);  

max_sum = 10000;
p1 = zeros(1,length_K);
% k=0
    k=0;
    i = 2:max_sum;
    sum_term = (-1).^i.*epsilon.^i./(factorial(i).*sqrt(i).*i);
    sum_term(find(isnan(sum_term)))=0;
    p1(1) = 1./(2*sqrt(2)*Q*(1+F)^2)...
        .*(sum(sum_term) - (1+F)*epsilon)...
        + 1;
% k=1
    k=1;
    i = 2:max_sum;
    sum_term = (-1).^i.*epsilon.^i./(factorial(i).*sqrt(i).*i);
    sum_term(find(isnan(sum_term)))=0;
    p1(2) = -1./(2*sqrt(2)*Q*(1+F)^2)...
    .*(sum(sum_term) - (1+F)*epsilon);
% k>=2
for k = 3:length_K
    real_k = k-1;
    i = real_k:max_sum;
    sum_term = (-1).^i.*epsilon.^i./(factorial(i-real_k).*sqrt(i).*i);
    sum_term(find(isnan(sum_term)))=0;
    p1(k) = (-1)^real_k/factorial(real_k)./(2*sqrt(2)*Q*(1+F)^2)...
    .*sum(sum_term);
end

% p1 = p1/sum(p1); % Normalize p1

N = 1E4; % max number of particles to tale into account. Should be higher than 1E3 since there are only 
% about 1.5 particles in the focus. Larger number of partices don't make the difference
pN = zeros(N,length_K);

% We compute convolutions (N) for values from k=1 to k=maxK
%This is equivalent to equation (20)
for i = 1:N
    if i > 2
        convp1pN = conv(p1,pN(i-1,:));
        pN(i,:) = convp1pN(1:length_K);   
    elseif i == 2
        pN(2,:) = p1; % First convolution
    elseif i == 1
        pN(1,1) = 1; % Probability of having 0 photons with 0 particles 
    end
end

% Compute equation (21)
vectorParticles = 0:N-1; % Note that we have to compute from 0 particles
PoisDiffusion = exp(-averageN).*averageN.^(vectorParticles)./factorial(vectorParticles);
PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0 (for high values of vectorParticles it take NaN values)

% Another option is to use a built in function but which is in the package
% of statistics: PoisDiffusion = poisspdf(vectorParticles, averageN);

% Compute equation (23)
PCH = zeros(1,length_K);
for k = 1:length_K
    PCH(k) = sum(pN(:,k).*PoisDiffusion');
end

% Normalize PCH acording to the sum to compare with the measurement
PCHnorm = PCH' / sum(PCH); % We have to put it vertical to match the shape of the vectors of gui_FCSfit