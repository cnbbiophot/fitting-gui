function PCHnorm=fitfcn_PCHsingleSpecies_normalized_single (param, tdata)

% PCH for ONE species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999
% does not work properly!
%
% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=averageN
% param(2)=epsilon

% agv, 15may2020

k = tdata;
maxK = length(k); % For the moment we dont take into account the k=0 value

averageN=param(1);    % Parámetros
epsilon=param(2);

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

PCHnorm = PCH' / sum(PCH); % We have to put it vertical to match the shape of the vectors of gui_FCSfit
