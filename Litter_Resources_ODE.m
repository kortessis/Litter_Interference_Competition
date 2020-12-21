function dXdt = Litter_Resources_ODE(t, X, w, a, beta, m, v, S, r, c)
% This is a funciton that specifies an ODE to pass into an ode solver in
% Matlab. We use ode45 in all cases.

% Vector to hold the growth rates in the ode
dXdt(:,1) = zeros(4,1);

% Plant sp1 growth equation
dXdt(1) = X(1).*(w(1)*a(1)*X(3)/(1+beta(1)*X(4)) - m(1)); 

% Plant sp2 growth equation
dXdt(2) = X(2).*(w(2)*a(2)*X(3)/(1+beta(2)*X(4)) - m(2));

% Resource Equation
dXdt(3) = r*(S - X(3)) - a(1)*X(3)*X(1)/(1 + beta(1)*X(4)) - a(2)*X(3)*X(2)/(1 + beta(2)*X(4));

% Litter Equation
dXdt(4) = v(1)*m(1)*X(1) + v(2)*m(2)*X(2) - c*X(4);
end
