function dXdt = Litter_ODE(t, X, r, beta, v, m, c)
% This is a function that creates the ode for passing into an ode solver in
% matlab.

% X is the species densities and litter
% r is the maximum growth rate when there is no litter, = waS - m
% beta is the sensitivity to litter
% v is the per-capita litter production rates
% c is the litter decomposition rate

n = length(r);
dXdt(:,1) = zeros(n+1,1);

dXdt(1:n) = X(1:n).*(r'./(1 + beta'*X(end)) - m');  % Plant growth equation
dXdt(end) = (v.*m)*X(1:n) - c*X(end);          % Litter Equation
end