function f=gaussbicurve (param, xdata)

% f=param(3)*gaussmf(xdata, [param(1), param(2)])+param(6)*gaussmf(xdata, [param(4), param(5)])+param(7);

% Returns the sum of two gausians, where 
% param(1) is sigma1, param(2) is the the centre, 
% param(3) is the amplitude for the first gaussian
% param(4) is sigma2, param(5) is the the centre, 
% param(6) is the amplitude for the second gaussian
% and param(7) is the baseline 



f=param(3)*exp(-0.5*((xdata-param(2))/param(1)).^2)+...
    param(6)*exp(-0.5*((xdata-param(5))/param(4)).^2)+param(7);