function f=gausscurve (param, xdata)

% Returns gausian, where param(1) is sigma, param(2) is the the centre, 
% and param(4) is the baseline and param(3) is the amplitude

%f = param(3)*gaussmf(xdata, [param(1), param(2)])+param(4);

f=param(3)*exp(-0.5*((xdata-param(2))/param(1)).^2)+param(4);