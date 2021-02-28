function F_finv = f_inv(P, eta1, eta2)

% F = f_inv(P, eta1, eta2) computes the inverse of the F_cdf with numerator 
% degrees of freedom eta1 and denominator degrees of freedom eta2 for the 
% corresponding probabilities in P. eta1 and eta2 parameters must contain real 
% positive values, and the values in P must lie on the interval [0 1].
%
% P corresponds to the probability that the f value lies OUT of the
% confidence interval (68% of probabilities that it is within the
% confidence interval implies, in our calculation P = 0.32)
% 
% It works ly and the table 4.4 (page 136) in Lakowicz can be
% reproduced. With the equivalence eta1=p (parameters fitted), 
% eta2=v (degrees of freedom) and P=P (probability that it is 
% out of the confidence interval).

% Task!!!! Implement an exception to return NaN when the function does not
% return a valid value (f_inv(0.025, 8, 169) returns a incorrect value
% while f_inv(0.025, 8, 160) returns the correct one)

% agv, 01Abr2019

lim_sup = 10000; % Lim sup of integration
num_int_points = 100000; % #Integration points
F_guess = 5; % This works fine

options = optimoptions('lsqnonlin');
F_finv = lsqnonlin(@err_curve, F_guess, [], [], options, P, eta1, eta2, lim_sup, num_int_points);

function err = err_curve (F, P, eta1, eta2, lim_sup, num_int_points)

f = linspace(F, lim_sup, num_int_points);

int_F = gamma((eta1+eta2)/2)/gamma(eta1/2)/gamma(eta2/2)*(eta1/eta2)^(eta1/2).*f.^((eta1-2)/2)./(1+f.*eta1/eta2).^((eta1+eta2)/2);
Pf = trapz(f,int_F);

err = (Pf - P);