function PCHnorm=fitfcn_PCHsingleSpecies_global (param, indParam, x)
           
% PCH for ONE species contained in the article "The Photon
% Counting Histogram in Fluorescence Fluctuation Spectroscopy". Y. Chen et
% al. 1999
%
% tdata es photon counts
%
% Definición de los parámetros actual: param(n)=Parameter Units
% Definición de los parámetros futura: param(n)=Parameter, Meaning, Units, Value, Lower bound, Upper bound
%
% param(1)=conc [M]
% param(2)=epsilon
% param(3)=z0 [m]
% param(4)=omega0 [m]

% agv, 15may2020
% agv, 22may2020

tdata=x(:,1);
% dsid=x(:,2); %We don't use dsid here because of integrations

number_sets = length(indParam(:,1));

k = tdata; % is a vector (1,2..., maxK, 1,2..., maxK, 1,2..., maxK , ...)
maxK = length(k)/number_sets; % For the moment we dont take into account the k=0 value

conc=param(indParam(:,1));    % Parámetros
epsilon=param(indParam(:,2));
z0=param(indParam(:,3));
omega0=param(indParam(:,4));

Na = 6.022E23;

VPSF = 4/3*pi*omega0.^2.*z0; % Volume of PSF in m^3

% Compute 1 particle probability with formula (16)
p1 = zeros(1,maxK*number_sets);

% we elliminate the for to introduce the dsid 
for ds = 1:number_sets %for the global fit
    PoisInt = @(x_fun) gammainc(abs(epsilon(ds)).*exp(-2*x_fun.^2),k(maxK*(ds-1)+1:maxK*ds)) .* gamma(k(maxK*(ds-1)+1:maxK*ds)); % Careful with the order of the input factors
    p1(maxK*(ds-1)+1:maxK*ds) = pi*omega0(ds).^2.*z0(ds)./(VPSF(ds) .* factorial(k(maxK*(ds-1)+1:maxK*ds))) .* integral(PoisInt,0,Inf,'ArrayValued',true); 
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

averageN = conc .* VPSF * 1E3 * Na; % Average number of particles
vectorParticles = 0:N-1;

for ds = 1:number_sets
    PoisDiffusion(ds,:) = exp(-averageN(ds)).*averageN(ds).^(vectorParticles)./...
        factorial(vectorParticles);
    PoisDiffusion(isnan(PoisDiffusion))=0; % Set all NaN to 0
end

PCH = zeros(1,maxK*number_sets);
PCHnorm = zeros(maxK*number_sets,1);
for ds = 1:number_sets
        PCH(maxK*(ds-1)+1:maxK*ds) = sum(pN(:,maxK*(ds-1)+1:maxK*ds).*PoisDiffusion(ds,:)')';
        PCHnorm(maxK*(ds-1)+1:maxK*ds,1) = PCH(maxK*(ds-1)+1:maxK*ds)/ sum(PCH(maxK*(ds-1)+1:maxK*ds));
end
